"""Calculations Browser tab: file browser with viewer, search, recalc, and more."""

import base64
import html as _html
import io
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import threading
import time
import zipfile
from collections import Counter
from itertools import permutations, product
from pathlib import Path

import ipywidgets as widgets
import numpy as np
from IPython.display import HTML, clear_output, display

from .constants import CALC_SEARCH_OPTIONS
from .energy_labels import format_xyz_comment_label
from .input_processing import (
    parse_inp_resources,
    parse_resource_settings,
    smiles_to_xyz,
    smiles_to_xyz_quick,
    contains_metal,
    is_smiles,
)
from .helpers import disable_spellcheck, save_neb_trajectory_csv, save_neb_trajectory_plot_png
from .molecule_viewer import (
    coord_to_xyz,
    parse_xyz_frames,
    DEFAULT_3DMOL_STYLE_JS,
    patch_viewer_mouse_controls_js,
)
from delfin.ensemble_nmr import CENSO_NMR_SOLVENT_CHOICES
from delfin.reporting.delfin_docx_report import _orca_plot_binary
from delfin.ssh_transfer_jobs import (
    create_transfer_job,
    ensure_jobs_dir,
    launch_transfer_job,
    list_transfer_jobs,
)
from delfin.user_settings import (
    get_settings_path,
    load_remote_archive_enabled,
    load_transfer_settings,
    normalize_ssh_transfer_settings,
)
from rdkit import Chem, RDLogger
from rdkit.Chem import rdDepictor, AllChem
from rdkit.Chem.Draw import MolToImage


ORCA_CPCM_TOP10_SOLVENTS = (
    "chloroform",
    "water",
    "acetonitrile",
    "dmso",
    "dmf",
    "methanol",
    "ethanol",
    "thf",
    "dichloromethane",
    "toluene",
)


def build_calc_nmr_input(coord_lines, *, pal: int, maxcore: int, solvent: str) -> str:
    """Build a standalone ORCA 1H NMR input with CPCM solvation."""
    solvent_name = str(solvent or "chloroform").strip() or "chloroform"
    coords_text = "\n".join(f"  {line.strip()}" for line in coord_lines if str(line).strip())
    return (
        f"!TPSS PCSSEG-1 AUTOAUX NMR CPCM({solvent_name})\n"
        f"\n"
        f"%pal\n"
        f"  nprocs {pal}\n"
        f"end\n"
        f"\n"
        f"%maxcore {maxcore}\n"
        f"\n"
        f"* xyz 0 1\n"
        f"{coords_text}\n"
        f"*\n"
        f"\n"
        f"%EPRNMR\n"
        f"   NUCLEI = ALL H {{SHIFT, SSALL}}\n"
        f"   TAU DOBSON\n"
        f"END\n"
    )


def create_tab(ctx):
    """Create the Calculations Browser tab.

    Returns ``(tab_widget, refs_dict)``.
    """
    # -- layout constants ---------------------------------------------------
    CALC_CONTENT_HEIGHT = 760
    CALC_MOL_SIZE = 450
    CALC_MOL_DYNAMIC_SCALE = 0.9725
    CALC_MOL_ZOOM = 0.9
    CALC_LEFT_DEFAULT = 375
    CALC_LEFT_MIN = 375
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
    CALC_XYZ_MAX_READ_BYTES = 50 * 1024 * 1024          # 50 MB – skip full read for huge trajectories
    CALC_CUBE_MAX_READ_BYTES = 100 * 1024 * 1024         # 100 MB – skip full read for huge cube files
    CALC_IMAGE_MAX_READ_BYTES = 20 * 1024 * 1024         # 20 MB – skip inline base64 for huge images
    CALC_LOG_XYZ_EXTRACT_MAX = 50 * 1024 * 1024          # 50 MB – skip XYZ extraction from huge logs
    # Trajectories with more frames than this use single-frame mode
    # (avoids embedding the full XYZ in a JS template literal via comm)
    CALC_XYZ_LARGE_TRAJ_FRAMES = 2000
    CALC_XYZ_PLAY_FPS_DEFAULT = 10
    CALC_XYZ_PLAY_FPS_MIN = 1
    CALC_XYZ_PLAY_FPS_MAX = 60
    CALC_VIEWER_PNG_SCALE = 6
    CALC_BROWSER_WORKFLOW_PAL = 4
    CALC_BROWSER_WORKFLOW_MAXCORE = 1000
    CALC_BROWSER_WORKFLOW_TIMELIMIT = '00:15:00'
    try:
        _calc_browser_root_hint = Path(
            os.environ.get('DELFIN_VOILA_ROOT_DIR') or str(Path.cwd())
        ).expanduser().resolve()
    except Exception:
        _calc_browser_root_hint = Path.cwd()
    CALC_BROWSER_UPLOAD_STAGING_REL = '.delfin_dashboard_uploads'
    CALC_BROWSER_UPLOAD_STAGING_DIR = _calc_browser_root_hint / CALC_BROWSER_UPLOAD_STAGING_REL
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
        'xyz_workflow_active': False,
        'xyz_workflow_mode': '',
        'xyz_workflow_source_path': None,
        'xyz_batch_path': '',
        'xyz_batch_marked_paths': [],
        'delete_current': False,
        'table_col_defs': [
            {'name': 'Value 1', 'type': 'text', 'pattern': '', 'occ': 'last'},
        ],
        'table_col_widgets': [],
        'table_csv_data': '',
        'table_panel_active': False,
        'file_is_preview': False,
        'file_preview_note': '',
        'file_chunk_start': 0,
        'file_chunk_end': 0,
        'chunk_dom_initialized': False,
        'xyz_frames': [],
        'xyz_current_frame': [0],
        'report_running': {},
        'traj_viewer_ready': False,
        'traj_playing': False,
        'traj_play_toggle_guard': False,
        'traj_play_stop_event': None,
        'rmsd_available': False,
        'rmsd_mode_active': False,
        'rmsd_saved_display': {},
        'rename_source_path': '',
        'duplicate_source_path': '',
        'ssh_transfer_running': False,
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
    upload_bridge_state = {
        'last_seq': 0,
        'files': {},
        'batches': {},
    }
    calc_scope_id = f'calc-scope-{abs(id(state))}'
    calc_resize_mol_fn = f'calcResizeMolViewer_{abs(id(state))}'
    calc_resize_pre_fn = f'calcResizePreselect3D_{abs(id(state))}'
    VIEWER_MOUSE_PATCH_JS = patch_viewer_mouse_controls_js('viewer', 'el')

    # -- widgets ------------------------------------------------------------
    calc_path_prefix = widgets.HTML(
        value='<b>📂 Path:</b>',
        layout=widgets.Layout(width='100%'),
    )
    calc_path_input = widgets.Text(
        value='/',
        continuous_update=False,
        layout=widgets.Layout(
            flex='1 1 0', min_width='0', width='1px', max_width='100%',
            height='24px', overflow_x='hidden', margin='0', padding='0',
        ),
    )
    calc_path_label = widgets.VBox(
        [calc_path_prefix, calc_path_input],
        layout=widgets.Layout(
            width='100%', overflow_x='hidden',
            align_items='stretch', gap='2px',
        ),
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
    calc_duplicate_btn = widgets.Button(
        description='Duplicate',
        layout=widgets.Layout(width='96px', height='26px'),
    )

    # Detect whether we are inside the Archive tab (calc_dir == archiv_dir)
    _is_archive_tab = ctx.calc_dir.resolve() == ctx.archive_dir.resolve()
    try:
        _remote_archive_enabled = load_remote_archive_enabled()
    except Exception:
        _remote_archive_enabled = False
    CALC_BROWSER_UPLOAD_SCOPE = 'archive' if _is_archive_tab else 'calculations'
    CALC_BROWSER_UPLOAD_STAGING_SCOPE_REL = (
        f'{CALC_BROWSER_UPLOAD_STAGING_REL}/{CALC_BROWSER_UPLOAD_SCOPE}'
    )
    CALC_BROWSER_UPLOAD_STAGING_SCOPE_DIR = (
        CALC_BROWSER_UPLOAD_STAGING_DIR / CALC_BROWSER_UPLOAD_SCOPE
    )

    # Move-to-Archive button (hidden when we are already inside the Archive)
    calc_move_archive_btn = widgets.Button(
        description='→ Archive', button_style='info',
        layout=widgets.Layout(
            width='90px', height='26px',
            display='none' if _is_archive_tab else 'inline-flex',
        ),
    )
    calc_back_to_calculations_btn = widgets.Button(
        description='← Calculations', button_style='info',
        layout=widgets.Layout(
            width='118px', height='26px',
            display='inline-flex' if _is_archive_tab else 'none',
        ),
    )
    calc_ssh_transfer_btn = widgets.Button(
        description='SSH Transfer',
        button_style='info',
        layout=widgets.Layout(
            width='110px', height='26px',
            display='inline-flex' if _remote_archive_enabled else 'none',
        ),
        disabled=not _remote_archive_enabled,
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

    # Explorer folder/file management
    calc_new_folder_input = widgets.Text(
        placeholder='New folder name',
        layout=widgets.Layout(flex='1 1 auto', min_width='120px', height='26px', display='none'),
    )
    calc_new_folder_btn = widgets.Button(
        description='+ Folder',
        button_style='success',
        layout=widgets.Layout(width='90px', height='26px', display='none'),
    )
    calc_rename_input = widgets.Text(
        placeholder='Rename selected/current folder',
        layout=widgets.Layout(flex='1 1 auto', min_width='120px', height='26px', display='none'),
    )
    calc_rename_btn = widgets.Button(
        description='Rename',
        layout=widgets.Layout(width='90px', height='26px', display='none'),
    )
    calc_move_target_input = widgets.Text(
        placeholder='Move target (relative or /from-root)',
        layout=widgets.Layout(flex='1 1 auto', min_width='120px', height='26px', display='none'),
    )
    calc_upload_target_input = widgets.Text(
        placeholder='Upload target (relative or /from-root)',
        layout=widgets.Layout(flex='1 1 auto', min_width='120px', height='26px', display='none'),
    )
    calc_action_source_input = widgets.Text(
        placeholder='Action source item (internal)',
        layout=widgets.Layout(flex='1 1 auto', min_width='120px', height='26px', display='none'),
    )
    calc_move_btn = widgets.Button(
        description='Move',
        layout=widgets.Layout(width='90px', height='26px', display='none'),
    )
    calc_hidden_upload = widgets.FileUpload(
        accept='',
        multiple=True,
        description='',
        layout=widgets.Layout(width='1px', height='1px', overflow='hidden'),
    )
    calc_upload_meta_input = widgets.Textarea(
        value='',
        layout=widgets.Layout(width='1px', height='1px', display='none'),
    )
    calc_upload_chunk_input = widgets.Textarea(
        value='',
        layout=widgets.Layout(width='1px', height='1px', display='none'),
    )
    calc_upload_seq_input = widgets.IntText(
        value=0,
        layout=widgets.Layout(width='1px', height='1px', display='none'),
    )
    calc_upload_ack_input = widgets.IntText(
        value=0,
        layout=widgets.Layout(width='1px', height='1px', display='none'),
    )
    calc_upload_trigger_btn = widgets.Button(
        description='',
        layout=widgets.Layout(width='1px', height='1px', display='none'),
    )
    calc_upload_ack_label = widgets.Label(
        value='0',
        layout=widgets.Layout(width='1px', height='1px', display='none'),
    )
    # Bridge widgets for JS → Python double-click and keyboard actions
    calc_dblclick_input = widgets.Text(
        value='',
        layout=widgets.Layout(width='1px', height='1px', display='none'),
    )
    calc_dblclick_input.add_class('calc-cmd-dblclick')
    calc_xyz_batch_dblclick_input = widgets.Text(
        value='',
        layout=widgets.Layout(width='1px', height='1px', display='none'),
    )
    calc_xyz_batch_dblclick_input.add_class('calc-cmd-xyz-batch-dblclick')
    calc_xyz_batch_toggle_input = widgets.Text(
        value='',
        layout=widgets.Layout(width='1px', height='1px', display='none'),
    )
    calc_xyz_batch_toggle_input.add_class('calc-cmd-xyz-batch-toggle')
    calc_keyboard_action_input = widgets.Text(
        value='',
        layout=widgets.Layout(width='1px', height='1px', display='none'),
    )
    calc_keyboard_action_input.add_class('calc-cmd-keyboard-action')
    calc_explorer_new_btn = widgets.Button(
        description='📁 New Folder',
        layout=widgets.Layout(width='110px', min_width='110px', height='26px'),
    )
    calc_explorer_rename_btn = widgets.Button(
        description='✏ Rename',
        layout=widgets.Layout(width='92px', height='26px'),
    )
    calc_rename_prompt_label = widgets.HTML('<b>Rename:</b>')
    calc_rename_prompt_input = widgets.Text(
        placeholder='New name',
        layout=widgets.Layout(flex='1 1 auto', min_width='120px', height='26px'),
    )
    calc_rename_prompt_ok_btn = widgets.Button(
        description='OK',
        button_style='primary',
        layout=widgets.Layout(width='56px', height='26px'),
    )
    calc_rename_prompt_cancel_btn = widgets.Button(
        description='Cancel',
        layout=widgets.Layout(width='78px', height='26px'),
    )
    calc_rename_prompt_row = widgets.HBox(
        [
            calc_rename_prompt_label,
            calc_rename_prompt_input,
            calc_rename_prompt_ok_btn,
            calc_rename_prompt_cancel_btn,
        ],
        layout=widgets.Layout(display='none', gap='6px', align_items='center'),
    )
    calc_duplicate_prompt_label = widgets.HTML('<b>Duplicate:</b>')
    calc_duplicate_prompt_input = widgets.Text(
        placeholder='Duplicate folder name',
        layout=widgets.Layout(flex='1 1 auto', min_width='120px', height='26px'),
    )
    calc_duplicate_prompt_yes_btn = widgets.Button(
        description='Yes',
        button_style='primary',
        layout=widgets.Layout(width='56px', height='26px'),
    )
    calc_duplicate_prompt_cancel_btn = widgets.Button(
        description='Cancel',
        layout=widgets.Layout(width='78px', height='26px'),
    )
    calc_duplicate_prompt_row = widgets.HBox(
        [
            calc_duplicate_prompt_label,
            calc_duplicate_prompt_input,
            calc_duplicate_prompt_yes_btn,
            calc_duplicate_prompt_cancel_btn,
        ],
        layout=widgets.Layout(display='none', gap='6px', align_items='center'),
    )
    calc_new_folder_prompt_label = widgets.HTML('<b>New folder name:</b>')
    calc_new_folder_prompt_input = widgets.Text(
        placeholder='Folder name',
        layout=widgets.Layout(flex='1 1 auto', min_width='120px', height='26px'),
    )
    calc_new_folder_prompt_ok_btn = widgets.Button(
        description='Create',
        button_style='primary',
        layout=widgets.Layout(width='70px', height='26px'),
    )
    calc_new_folder_prompt_cancel_btn = widgets.Button(
        description='Cancel',
        layout=widgets.Layout(width='78px', height='26px'),
    )
    calc_new_folder_prompt_row = widgets.HBox(
        [
            calc_new_folder_prompt_label,
            calc_new_folder_prompt_input,
            calc_new_folder_prompt_ok_btn,
            calc_new_folder_prompt_cancel_btn,
        ],
        layout=widgets.Layout(display='none', gap='6px', align_items='center'),
    )
    calc_transfer_jobs_btn = widgets.Button(
        description='Running Transfers',
        layout=widgets.Layout(
            width='140px', height='26px',
            display='inline-flex' if _remote_archive_enabled else 'none',
        ),
        disabled=not _remote_archive_enabled,
    )
    calc_transfer_jobs_refresh_btn = widgets.Button(
        description='Refresh Jobs',
        layout=widgets.Layout(width='106px', height='26px'),
    )
    calc_transfer_jobs_html = widgets.HTML(
        value='',
        layout=widgets.Layout(
            width='100%',
            overflow_x='hidden',
            overflow_y='auto',
            flex='1 1 0',
            min_height='0',
        ),
    )
    calc_transfer_jobs_panel = widgets.VBox(
        [
            widgets.HBox(
                [
                    widgets.HTML('<b>Running Transfers</b>'),
                    calc_transfer_jobs_refresh_btn,
                ],
                layout=widgets.Layout(gap='6px', align_items='center', flex_flow='row wrap'),
            ),
            calc_transfer_jobs_html,
        ],
        layout=widgets.Layout(
            display='none',
            gap='6px',
            width='100%',
            flex='1 1 0',
            min_height='0',
            overflow='hidden',
        ),
    )
    calc_ops_status = widgets.HTML(
        value='',
        layout=widgets.Layout(width='100%', overflow_x='hidden'),
    )
    calc_path_label.add_class('calc-path-label')
    calc_path_input.add_class('calc-cmd-path-input')
    calc_path_input.add_class('delfin-nospell')
    calc_refresh_btn.add_class('calc-cmd-refresh-btn')
    calc_new_folder_input.add_class('calc-cmd-new-folder-input')
    calc_new_folder_btn.add_class('calc-cmd-new-folder-btn')
    calc_explorer_new_btn.add_class('calc-cmd-explorer-new-btn')
    calc_rename_input.add_class('calc-cmd-rename-input')
    calc_rename_btn.add_class('calc-cmd-rename-btn')
    calc_move_target_input.add_class('calc-cmd-move-target')
    calc_upload_target_input.add_class('calc-cmd-upload-target')
    calc_action_source_input.add_class('calc-cmd-action-source')
    calc_move_btn.add_class('calc-cmd-move-btn')
    calc_hidden_upload.add_class('calc-hidden-upload')
    calc_ops_status.add_class('calc-ops-status')
    calc_upload_meta_input.add_class('calc-upload-meta')
    calc_upload_chunk_input.add_class('calc-upload-chunk')
    calc_upload_seq_input.add_class('calc-upload-seq')
    calc_upload_ack_input.add_class('calc-upload-ack')
    calc_upload_trigger_btn.add_class('calc-upload-trigger-btn')
    calc_upload_ack_label.add_class('calc-upload-ack-label')

    # ---- Table extraction panel (Archive tab only) -------------------------
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

    def _normalize_table_number(token):
        value = str(token).strip().replace('D', 'E').replace('d', 'e')
        # Accept decimal comma from some text exports.
        if ',' in value and '.' not in value:
            value = value.replace(',', '.')
        return value

    calc_table_btn = widgets.Button(
        description='📊 Table',
        layout=widgets.Layout(
            width='84px', height='26px',
            display='inline-flex' if _is_archive_tab else 'none',
        ),
    )
    calc_table_file_input = widgets.Text(
        placeholder='e.g. orca.out or result.json',
        layout=widgets.Layout(flex='1 1 auto', min_width='120px', height='26px'),
    )
    calc_table_scope_dd = widgets.Dropdown(
        options=[('All folders', 'all'), ('Selection only', 'selected')],
        value='all',
        layout=widgets.Layout(width='130px', height='26px'),
    )
    calc_table_recursive_cb = widgets.ToggleButton(
        value=False, description='Recurse',
        layout=widgets.Layout(width='86px', height='26px'),
    )
    calc_table_decimal_comma_btn = widgets.ToggleButton(
        value=False, description='Dot to Comma',
        tooltip='Convert decimal point to decimal comma in table/CSV values',
        layout=widgets.Layout(width='126px', height='26px'),
    )
    calc_table_preset_name = widgets.Text(
        placeholder='extract_table.delfin-table.json',
        layout=widgets.Layout(flex='1 1 auto', min_width='140px', height='26px'),
    )
    calc_table_preset_name.add_class('calc-table-preset-name')
    calc_table_preset_save_btn = widgets.Button(
        description='Save file',
        layout=widgets.Layout(width='88px', height='26px'),
    )
    calc_table_preset_status = widgets.HTML(
        value="<span style='color:#616161;'>Click a <code>.delfin-table.json</code> file in the explorer to load it.</span>"
    )
    calc_table_add_col_btn = widgets.Button(
        description='+ Column',
        layout=widgets.Layout(width='80px', height='26px'),
    )
    calc_table_run_btn = widgets.Button(
        description='▶ Run', button_style='primary',
        layout=widgets.Layout(width='68px', height='26px'),
    )
    calc_table_close_btn = widgets.Button(
        description='✕',
        layout=widgets.Layout(width='32px', height='26px'),
    )
    calc_table_csv_btn = widgets.Button(
        description='⬇ CSV',
        layout=widgets.Layout(width='68px', height='26px', display='none'),
    )
    calc_table_status = widgets.HTML(value='')
    calc_table_output = widgets.HTML(
        value='',
        layout=widgets.Layout(
            width='100%', overflow_x='auto', overflow_y='auto',
            flex='1 1 0', min_height='0',
        ),
    )
    calc_table_cols_box = widgets.VBox(
        [], layout=widgets.Layout(width='100%', gap='4px'),
    )
    calc_table_file_row = widgets.HBox([
        widgets.HTML('<span style="white-space:nowrap"><b>File:</b></span>'),
        calc_table_file_input,
        calc_table_scope_dd,
        calc_table_recursive_cb,
        calc_table_decimal_comma_btn,
    ], layout=widgets.Layout(gap='6px', align_items='center', width='100%'))
    calc_table_file_row.add_class('calc-table-top-row')

    calc_table_preset_row = widgets.HBox([
        widgets.HTML('<span style="white-space:nowrap"><b>Preset file:</b></span>'),
        calc_table_preset_name,
        calc_table_preset_save_btn,
    ], layout=widgets.Layout(gap='6px', align_items='center', width='100%'))
    calc_table_preset_row.add_class('calc-table-top-row')

    calc_table_panel = widgets.VBox([
        widgets.HBox([
            widgets.HTML('<b>📊 Extract Table</b>'),
            calc_table_close_btn,
        ], layout=widgets.Layout(
            justify_content='space-between', align_items='center', width='100%',
        )),
        calc_table_file_row,
        calc_table_preset_row,
        calc_table_preset_status,
        calc_table_cols_box,
        widgets.HBox(
            [calc_table_add_col_btn, calc_table_run_btn, calc_table_csv_btn],
            layout=widgets.Layout(gap='6px', align_items='center'),
        ),
        calc_table_status,
        calc_table_output,
    ], layout=widgets.Layout(
        display='none', width='100%', padding='8px', gap='6px',
        border='1px solid #e0e0e0', border_radius='4px', overflow_x='hidden',
        flex='1 1 0', min_height='0',
    ))

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

    # File list (always multi-select; JS handles click semantics)
    calc_file_list = widgets.SelectMultiple(
        options=[], rows=22,
        layout=widgets.Layout(
            width='100%', flex='1 1 0', min_height='0', margin='-4px 0 0 0'
        ),
    )
    calc_file_list.add_class('calc-file-list')

    # Content area
    calc_content_area = widgets.HTML(
        value='',
        layout=widgets.Layout(
            width='100%', display='block', overflow_x='hidden',
            flex='1 1 0', min_height='0',
        ),
    )

    # Molecule viewer
    calc_view_png_btn = widgets.Button(
        description='PNG',
        icon='',
        button_style='warning',
        layout=widgets.Layout(width='54px', min_width='54px', height='30px', display='none'),
        tooltip='Download a high-resolution PNG from the current 3D view',
    )
    calc_view_png_btn.add_class('calc-png-btn')
    calc_mol_label = widgets.HTML(
        "<div style='height:26px; line-height:26px; margin:0;'>"
        "<b>🔬 Molecule Preview:</b></div>"
    )
    calc_mol_header = widgets.HBox(
        [calc_mol_label, calc_view_png_btn],
        layout=widgets.Layout(
            width='100%',
            align_items='flex-start',
            justify_content='space-between',
            margin='0 0 8px 0',
        ),
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
        value='', layout=widgets.Layout(margin='0 0 6px 0'),
    )
    calc_xyz_frame_input = widgets.BoundedIntText(
        value=1, min=1, max=1, step=1,
        layout=widgets.Layout(width='72px', min_width='72px', height='28px'),
    )
    calc_xyz_frame_input.add_class('calc-xyz-frame-input')
    calc_xyz_frame_total = widgets.HTML(
        value='<b>/ 1</b>', layout=widgets.Layout(width='62px', min_width='62px'),
    )
    calc_xyz_loop_checkbox = widgets.ToggleButton(
        value=True,
        description='Loop',
        button_style='info',
        layout=widgets.Layout(width='82px', height='32px'),
    )
    calc_xyz_fps_input = widgets.BoundedIntText(
        value=CALC_XYZ_PLAY_FPS_DEFAULT,
        min=CALC_XYZ_PLAY_FPS_MIN,
        max=CALC_XYZ_PLAY_FPS_MAX,
        step=1,
        layout=widgets.Layout(width='72px', height='28px'),
    )
    calc_xyz_play_btn = widgets.ToggleButton(
        value=False,
        description='Play',
        icon='play',
        button_style='success',
        layout=widgets.Layout(width='86px', height='32px'),
    )
    calc_xyz_play_btn.add_class('calc-xyz-play-btn')
    calc_xyz_copy_btn = widgets.Button(
        description='📋 Copy Coordinates', button_style='success',
        layout=widgets.Layout(width='176px', min_width='176px', height='32px'),
    )
    calc_xyz_png_btn = widgets.Button(
        description='🖼 PNG', button_style='warning',
        layout=widgets.Layout(width='112px', min_width='112px', height='32px'),
        tooltip='Download a high-resolution PNG from the current XYZ view',
    )
    calc_xyz_frame_inline = widgets.HBox(
        [widgets.HTML('<b>Frame:</b>'), calc_xyz_frame_input, calc_xyz_frame_total],
        layout=widgets.Layout(gap='10px', align_items='center', min_width='170px', flex='0 0 auto'),
    )
    calc_xyz_png_row = widgets.HBox(
        [calc_xyz_png_btn],
        layout=widgets.Layout(
            display='none',
            width='100%',
            justify_content='flex-end',
            align_items='center',
            margin='0 0 6px 0',
        ),
    )
    calc_xyz_controls = widgets.HBox(
        [
            calc_xyz_frame_inline,
            calc_xyz_copy_btn,
        ],
        layout=widgets.Layout(
            display='none',
            gap='12px',
            margin='0 0 6px 0',
            align_items='center',
            justify_content='space-between',
            flex_flow='row nowrap',
            width='100%',
        ),
    )
    # Coord (Turbomole) copy button
    calc_coord_copy_btn = widgets.Button(
        description='📋 Copy Coordinates', button_style='success',
        layout=widgets.Layout(width='160px', height='32px'),
    )
    calc_coord_png_btn = widgets.Button(
        description='🖼 PNG', button_style='warning',
        layout=widgets.Layout(width='108px', height='32px'),
        tooltip='Download a high-resolution PNG from the current 3D view',
    )
    calc_coord_png_row = widgets.HBox(
        [calc_coord_png_btn],
        layout=widgets.Layout(
            display='none',
            width='100%',
            justify_content='flex-end',
            align_items='center',
            margin='0 0 6px 0',
        ),
    )
    calc_coord_controls = widgets.HBox(
        [calc_coord_copy_btn],
        layout=widgets.Layout(
            display='none', gap='10px', margin='8px 0',
            align_items='center', justify_content='flex-end', width='100%',
        ),
    )
    calc_xyz_tray_controls = widgets.VBox(
        [
            calc_xyz_frame_label,
            calc_xyz_png_row,
            calc_xyz_controls,
            calc_coord_png_row,
            calc_coord_controls,
            widgets.HBox(
                [
                    calc_xyz_loop_checkbox,
                    widgets.HBox(
                        [widgets.HTML('<b>FPS:</b>'), calc_xyz_fps_input],
                        layout=widgets.Layout(gap='6px', align_items='center', flex='0 0 auto'),
                    ),
                    calc_xyz_play_btn,
                ],
                layout=widgets.Layout(
                    gap='12px', align_items='center', width='100%', justify_content='space-between',
                ),
            ),
        ],
        layout=widgets.Layout(
            display='none',
            gap='14px',
            align_items='stretch',
            width='360px',
            min_width='340px',
            max_width='420px',
            margin='0',
        ),
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
    calc_xyz_workflow_info = widgets.HTML(
        value='', layout=widgets.Layout(width='100%', overflow_x='hidden'),
    )
    calc_xyz_workflow_pal = widgets.BoundedIntText(
        value=CALC_BROWSER_WORKFLOW_PAL,
        min=1,
        max=999,
        description='PAL',
        layout=widgets.Layout(width='150px', height='26px'),
    )
    calc_xyz_workflow_time = widgets.Text(
        value=CALC_BROWSER_WORKFLOW_TIMELIMIT,
        description='JobTime',
        layout=widgets.Layout(width='220px', height='26px'),
    )
    calc_xyz_workflow_bfw = widgets.ToggleButton(
        value=False,
        description='BFW Off',
        tooltip='Aktiviere -BFW fuer tadf_xtb/std2.',
        button_style='',
        layout=widgets.Layout(width='110px', min_width='110px', height='26px', display='none'),
    )
    calc_submit_xyz_workflow_btn = widgets.Button(
        description='Submit', button_style='success',
        layout=widgets.Layout(width='90px', min_width='90px', height='26px'), disabled=True,
    )
    calc_xyz_workflow_status = widgets.HTML(
        value='', layout=widgets.Layout(width='100%', overflow_x='hidden'),
    )

    def _calc_update_bfw_toggle(*_args):
        enabled = bool(calc_xyz_workflow_bfw.value)
        calc_xyz_workflow_bfw.description = 'BFW On' if enabled else 'BFW Off'
        calc_xyz_workflow_bfw.button_style = 'success' if enabled else ''

    calc_xyz_workflow_bfw.observe(_calc_update_bfw_toggle, names='value')
    _calc_update_bfw_toggle()
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
        value='', placeholder='STAGE=INDEX[,STAGE=INDEX,...]',
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

    # Build Batch from XYZ widgets
    calc_xyz_batch_filename = widgets.Text(
        value='',
        placeholder='batch_name',
        layout=widgets.Layout(width='220px', min_width='160px', height='26px'),
    )
    calc_xyz_batch_refresh_btn = widgets.Button(
        description='Refresh',
        layout=widgets.Layout(width='88px', min_width='88px', height='26px'),
    )
    calc_xyz_batch_up_btn = widgets.Button(
        description='Up',
        layout=widgets.Layout(width='70px', min_width='70px', height='26px'),
    )
    calc_xyz_batch_root_btn = widgets.Button(
        description='Calc',
        layout=widgets.Layout(width='70px', min_width='70px', height='26px'),
    )
    calc_xyz_batch_mark_btn = widgets.Button(
        description='Select',
        button_style='warning',
        layout=widgets.Layout(width='80px', min_width='80px', height='26px'),
    )
    calc_xyz_batch_select_all_btn = widgets.Button(
        description='Select All',
        layout=widgets.Layout(width='95px', min_width='95px', height='26px'),
    )
    calc_xyz_batch_build_btn = widgets.Button(
        description='Build Batch TXT', button_style='primary',
        layout=widgets.Layout(width='130px', min_width='130px', height='26px'),
    )
    calc_xyz_batch_copy_btn = widgets.Button(
        description='Copy Batch TXT', button_style='success',
        layout=widgets.Layout(width='130px', min_width='130px', height='26px'),
    )
    calc_xyz_batch_select = widgets.SelectMultiple(
        options=[], rows=10,
        layout=widgets.Layout(width='100%', min_height='170px', max_height='300px'),
    )
    calc_xyz_batch_select.add_class('delfin-multi-select')
    calc_xyz_batch_select.add_class('calc-xyz-batch-select')
    calc_xyz_batch_root_info = widgets.HTML(
        value='', layout=widgets.Layout(width='100%', overflow_x='hidden'),
    )
    calc_xyz_batch_status = widgets.HTML(
        value='', layout=widgets.Layout(width='100%', overflow_x='hidden'),
    )
    calc_xyz_batch_panel = widgets.VBox(
        [
            widgets.HTML('<b>Build Batch from XYZ</b>'),
            widgets.HTML(
                '<span style="color:#555;">Click left box to set/remove Haken. Double-click folders to enter them.</span>'
            ),
            widgets.HBox(
                [
                    calc_xyz_batch_up_btn,
                    calc_xyz_batch_root_btn,
                    calc_xyz_batch_mark_btn,
                    calc_xyz_batch_select_all_btn,
                    calc_xyz_batch_root_info,
                ],
                layout=widgets.Layout(
                    gap='8px', align_items='center', flex_flow='row wrap', width='100%',
                ),
            ),
            widgets.HBox(
                [
                    widgets.HTML('<b>Batch_name:</b>'),
                    calc_xyz_batch_filename,
                    calc_xyz_batch_refresh_btn,
                    calc_xyz_batch_build_btn,
                    calc_xyz_batch_copy_btn,
                ],
                layout=widgets.Layout(
                    gap='8px', align_items='center', flex_flow='row wrap', width='100%',
                ),
            ),
            calc_xyz_batch_select,
            calc_xyz_batch_status,
        ],
        layout=widgets.Layout(display='none', margin='8px 0 8px 0', width='100%'),
    )

    # Print Mode widgets (.out -> .hess via orca_pltvib)
    calc_print_mode_input = widgets.Text(
        value='',
        placeholder='e.g. 6 7 8 or 6,7,8',
        layout=widgets.Layout(width='220px', min_width='180px', height='26px'),
    )
    calc_print_mode_plot_btn = widgets.Button(
        description='Plot',
        button_style='primary',
        layout=widgets.Layout(width='88px', min_width='88px', height='26px'),
    )
    calc_print_mode_file_info = widgets.HTML(
        value='',
        layout=widgets.Layout(width='100%', overflow_x='hidden'),
    )
    calc_print_mode_status = widgets.HTML(
        value='',
        layout=widgets.Layout(width='100%', overflow_x='hidden'),
    )
    calc_print_mode_panel = widgets.VBox(
        [
            widgets.HTML('<b>Print Mode</b>'),
            widgets.HTML(
                '<span style="color:#555;">Enter modes to print.</span>'
            ),
            calc_print_mode_file_info,
            widgets.HBox(
                [widgets.HTML('<b>Modes:</b>'), calc_print_mode_input, calc_print_mode_plot_btn],
                layout=widgets.Layout(
                    gap='8px', align_items='center', flex_flow='row wrap', width='100%',
                ),
            ),
            calc_print_mode_status,
        ],
        layout=widgets.Layout(display='none', margin='8px 0 8px 0', width='100%'),
    )

    # MO Plot widgets (.out/.gbw -> .cube via orca_plot)
    calc_mo_plot_input = widgets.Text(
        value='',
        placeholder='e.g. 45 or 45,46,47',
        layout=widgets.Layout(width='180px', min_width='160px', height='26px'),
    )
    calc_mo_plot_alpha_cb = widgets.ToggleButton(
        value=True,
        description='Alpha',
        tooltip='Plot alpha orbitals',
        layout=widgets.Layout(width='86px', min_width='86px', height='26px'),
    )
    calc_mo_plot_beta_cb = widgets.ToggleButton(
        value=True,
        description='Beta',
        tooltip='Plot beta orbitals',
        layout=widgets.Layout(width='86px', min_width='86px', height='26px'),
    )
    calc_mo_plot_spin_box = widgets.HBox(
        [calc_mo_plot_alpha_cb, calc_mo_plot_beta_cb],
        layout=widgets.Layout(display='none', gap='6px'),
    )
    calc_mo_plot_btn = widgets.Button(
        description='Plot MO',
        button_style='primary',
        layout=widgets.Layout(width='100px', min_width='100px', height='26px'),
    )
    calc_mo_plot_file_info = widgets.HTML(
        value='',
        layout=widgets.Layout(width='100%', overflow_x='hidden'),
    )
    calc_mo_plot_status = widgets.HTML(
        value='',
        layout=widgets.Layout(width='100%', overflow_x='hidden'),
    )
    calc_mo_plot_panel = widgets.VBox(
        [
            widgets.HTML('<b>MO Plot</b>'),
            widgets.HTML(
                '<span style="color:#555;">Generate orbital cube files via orca_plot and display in 3D viewer.</span>'
            ),
            calc_mo_plot_file_info,
            widgets.HBox(
                [
                    widgets.HTML('<b>MO #:</b>'),
                    calc_mo_plot_input,
                    calc_mo_plot_spin_box,
                    calc_mo_plot_btn,
                ],
                layout=widgets.Layout(
                    gap='8px', align_items='center', flex_flow='row wrap', width='100%',
                ),
            ),
            calc_mo_plot_status,
        ],
        layout=widgets.Layout(display='none', margin='8px 0 8px 0', width='100%'),
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
    calc_mol_view_wrap = widgets.Box(
        [calc_mol_viewer],
        layout=widgets.Layout(flex='0 0 auto', min_width='0', width='auto'),
    )
    calc_mol_view_row = widgets.HBox(
        [calc_mol_view_wrap, calc_xyz_tray_controls],
        layout=widgets.Layout(
            width='100%',
            gap='12px',
            align_items='flex-start',
            justify_content='flex-start',
            flex_flow='row nowrap',
        ),
    )
    calc_mol_view_row.add_class('calc-mol-view-row')
    calc_mol_view_wrap.add_class('calc-mol-view-wrap')
    calc_xyz_tray_controls.add_class('calc-xyz-tray-controls')
    calc_mol_container = widgets.VBox(
        [
            calc_mol_header,
            calc_rmsd_controls,
            calc_mol_view_row,
        ],
        layout=widgets.Layout(display='none', margin='0 0 10px 0', width='100%', align_items='stretch'),
    )

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
    calc_xyz_workflow_toolbar = widgets.VBox([
        calc_xyz_workflow_info,
        widgets.HBox(
            [calc_xyz_workflow_pal, calc_xyz_workflow_time, calc_xyz_workflow_bfw, calc_submit_xyz_workflow_btn],
            layout=widgets.Layout(
                width='100%', overflow_x='hidden', gap='8px',
                flex_flow='row wrap', align_items='center',
            ),
        ),
        calc_xyz_workflow_status,
    ], layout=widgets.Layout(
        margin='5px 0', width='100%', overflow_x='hidden', gap='6px',
        display='none',
    ))
    # --- Calc NMR panel widgets ---
    calc_nmr_pal = widgets.BoundedIntText(
        value=12, min=1, max=999, description='PAL',
        layout=widgets.Layout(width='130px', height='26px'),
    )
    calc_nmr_maxcore = widgets.BoundedIntText(
        value=3000, min=100, max=99999, description='MaxCore',
        layout=widgets.Layout(width='160px', height='26px'),
    )
    calc_nmr_time = widgets.Text(
        value='03:00:00', description='JobTime',
        layout=widgets.Layout(width='200px', height='26px'),
    )
    calc_nmr_solvent = widgets.Dropdown(
        options=list(ORCA_CPCM_TOP10_SOLVENTS),
        value='chloroform',
        description='Solvent',
        layout=widgets.Layout(width='220px', height='26px'),
    )
    calc_nmr_submit_btn = widgets.Button(
        description='Submit', button_style='primary',
        layout=widgets.Layout(width='100px', height='26px'),
    )
    calc_nmr_status = widgets.HTML(value='')
    calc_nmr_panel = widgets.VBox([
        widgets.HTML('<b>Calc NMR</b>'),
        widgets.HTML(
            '<span style="color:#555;">Generate ORCA 1H NMR input and submit SLURM job. '
            'Creates an NMR/ subfolder with the .inp file.</span>'
        ),
        widgets.HBox(
            [calc_nmr_pal, calc_nmr_maxcore, calc_nmr_time, calc_nmr_solvent, calc_nmr_submit_btn],
            layout=widgets.Layout(
                width='100%', overflow_x='hidden', gap='8px',
                flex_flow='row wrap', align_items='center',
            ),
        ),
        calc_nmr_status,
    ], layout=widgets.Layout(
        display='none', margin='8px 0 8px 0', width='100%',
    ))
    # --- Calc CENSO/ANMR panel widgets ---
    calc_censo_nmr_pal = widgets.BoundedIntText(
        value=12, min=1, max=999, description='PAL',
        layout=widgets.Layout(width='130px', height='26px'),
    )
    calc_censo_nmr_maxcore = widgets.BoundedIntText(
        value=3000, min=100, max=99999, description='MaxCore',
        layout=widgets.Layout(width='160px', height='26px'),
    )
    calc_censo_nmr_time = widgets.Text(
        value='12:00:00', description='JobTime',
        layout=widgets.Layout(width='200px', height='26px'),
    )
    calc_censo_nmr_solvent = widgets.Dropdown(
        options=list(CENSO_NMR_SOLVENT_CHOICES),
        value='chcl3',
        description='Solvent',
        layout=widgets.Layout(width='220px', height='26px'),
    )
    calc_censo_nmr_charge = widgets.BoundedIntText(
        value=0, min=-99, max=99, description='Charge',
        layout=widgets.Layout(width='150px', height='26px'),
    )
    calc_censo_nmr_multiplicity = widgets.BoundedIntText(
        value=1, min=1, max=20, description='Mult',
        layout=widgets.Layout(width='140px', height='26px'),
    )
    calc_censo_nmr_mhz = widgets.BoundedIntText(
        value=400, min=50, max=2000, description='MHz',
        layout=widgets.Layout(width='140px', height='26px'),
    )
    calc_censo_nmr_submit_btn = widgets.Button(
        description='Submit', button_style='primary',
        layout=widgets.Layout(width='100px', height='26px'),
    )
    calc_censo_nmr_status = widgets.HTML(value='')
    calc_censo_nmr_panel = widgets.VBox([
        widgets.HTML('<b>Calc CENSO/ANMR</b>'),
        widgets.HTML(
            '<span style="color:#555;">Run an ensemble 1H NMR workflow via '
            'CREST + CENSO + ANMR in a new subfolder. Missing CENSO/ANMR helper '
            'tools are auto-installed on demand; crest, xtb, and ORCA must already exist.</span>'
        ),
        widgets.HBox(
            [
                calc_censo_nmr_pal,
                calc_censo_nmr_maxcore,
                calc_censo_nmr_time,
                calc_censo_nmr_solvent,
                calc_censo_nmr_charge,
                calc_censo_nmr_multiplicity,
                calc_censo_nmr_mhz,
                calc_censo_nmr_submit_btn,
            ],
            layout=widgets.Layout(
                width='100%', overflow_x='hidden', gap='8px',
                flex_flow='row wrap', align_items='center',
            ),
        ),
        calc_censo_nmr_status,
    ], layout=widgets.Layout(
        display='none', margin='8px 0 8px 0', width='100%',
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
        label_text = str(label or '').replace('\xa0', ' ').strip()
        if not label_text or label_text.startswith('('):
            return ''
        while True:
            parts = re.split(r'\s+', label_text, maxsplit=1)
            if len(parts) != 2:
                return label_text
            head, tail = parts[0].strip(), parts[1].strip()
            if head in {'📂', '🔬', '✅', '⬜', '☑', '☐', '✔', '□', '📁', '📄'} or not any(ch.isalnum() for ch in head):
                label_text = tail
                continue
            return label_text

    def _calc_path_for_label(label):
        name = _calc_label_to_name(label)
        if not name:
            return None
        return (
            (_calc_dir() / state['current_path'] / name)
            if state['current_path']
            else (_calc_dir() / name)
        )

    def _calc_resolve_active_label(preferred_label=''):
        candidates = []
        selected_labels = _calc_selected_labels()
        if len(selected_labels) == 1:
            candidates.append(selected_labels[0])
        if preferred_label and preferred_label not in candidates:
            candidates.append(preferred_label)

        for label in candidates:
            path = _calc_path_for_label(label)
            if path is not None and path.exists():
                return label, path

        for label in candidates:
            path = _calc_path_for_label(label)
            if path is not None:
                return label, path

        return '', None

    def _calc_get_selected_path():
        labels = _calc_selected_labels()
        if not labels:
            return None
        return _calc_path_for_label(labels[0])

    def _calc_set_ops_status(message, color='#555'):
        calc_ops_status.value = f'<span style="color:{color};">{message}</span>'

    def _calc_current_dir():
        return _calc_dir() / state['current_path'] if state['current_path'] else _calc_dir()

    def _calc_state_rel_path(path):
        rel = Path(path).resolve().relative_to(_calc_dir().resolve()).as_posix()
        return '' if rel in {'', '.'} else rel

    def _calc_label_for_name(name):
        target_name = str(name or '').strip()
        if not target_name:
            return ''
        for label in list(state.get('all_items') or []):
            if _calc_label_to_name(label) == target_name:
                return label
        return ''

    def _calc_safe_job_token(text, fallback):
        cleaned = re.sub(r'[^A-Za-z0-9._-]+', '_', str(text or '').strip()).strip('._-')
        return cleaned or fallback

    def _calc_count_xyz_frames(path, stop_after=2):
        frames = 0
        try:
            with Path(path).open('r', encoding='utf-8', errors='ignore') as handle:
                while frames < max(1, int(stop_after)):
                    atom_line = None
                    for raw_line in handle:
                        stripped = raw_line.strip()
                        if stripped:
                            atom_line = stripped
                            break
                    if atom_line is None:
                        break
                    try:
                        atom_count = int(atom_line)
                    except ValueError:
                        return 0
                    if atom_count <= 0:
                        return 0
                    if next(handle, None) is None:
                        return 0
                    for _ in range(atom_count):
                        if next(handle, None) is None:
                            return 0
                    frames += 1
        except OSError:
            return 0
        return frames

    def _calc_is_single_structure_xyz(path):
        candidate = Path(path)
        return (
            candidate.is_file()
            and candidate.suffix.lower() == '.xyz'
            and _calc_count_xyz_frames(candidate, stop_after=2) == 1
        )

    def _calc_next_available_dir(parent_dir, folder_name):
        candidate = Path(parent_dir) / folder_name
        if not candidate.exists():
            return candidate
        index = 2
        while True:
            numbered = Path(parent_dir) / f'{folder_name}_{index}'
            if not numbered.exists():
                return numbered
            index += 1

    def _calc_next_available_path(parent_dir, filename):
        base_name = str(filename or '').replace('\\', '/').split('/')[-1]
        base_name = _calc_clean_item_name(base_name)
        candidate = Path(parent_dir) / base_name
        if not candidate.exists():
            return candidate
        suffix = ''.join(Path(base_name).suffixes)
        stem = base_name[:-len(suffix)] if suffix else base_name
        index = 2
        while True:
            numbered_name = f'{stem}_{index}{suffix}'
            numbered = Path(parent_dir) / numbered_name
            if not numbered.exists():
                return numbered
            index += 1

    def _calc_reset_options_dropdown():
        try:
            calc_options_dropdown.unobserve(calc_on_options_change, names='value')
        except Exception:
            pass
        calc_options_dropdown.value = '(Options)'
        calc_options_dropdown.observe(calc_on_options_change, names='value')

    def _calc_prepare_xyz_browser_workflow(mode):
        selected_path = _calc_selected_item_path()
        if selected_path is None or selected_path.suffix.lower() != '.xyz':
            calc_reset_xyz_workflow_state()
            calc_update_view()
            _calc_set_ops_status('Select a single-structure .xyz file first.', '#d32f2f')
            return
        if not _calc_is_single_structure_xyz(selected_path):
            calc_reset_xyz_workflow_state()
            calc_update_view()
            _calc_set_ops_status('This action is only available for .xyz files with exactly one structure.', '#d32f2f')
            return

        if mode not in ('hyperpol_xtb', 'tadf_xtb'):
            calc_reset_xyz_workflow_state()
            calc_update_view()
            _calc_set_ops_status(f'Unknown workflow mode: {_html.escape(str(mode))}', '#d32f2f')
            return

        state['xyz_workflow_active'] = True
        state['xyz_workflow_mode'] = mode
        state['xyz_workflow_source_path'] = selected_path
        calc_xyz_workflow_pal.value = CALC_BROWSER_WORKFLOW_PAL
        calc_xyz_workflow_time.value = CALC_BROWSER_WORKFLOW_TIMELIMIT
        calc_xyz_workflow_bfw.layout.display = 'inline-flex'
        calc_submit_xyz_workflow_btn.disabled = False
        calc_xyz_workflow_info.value = (
            f'<b>{_html.escape(mode)}</b> for '
            f'<code>{_html.escape(selected_path.name)}</code>'
        )
        calc_xyz_workflow_status.value = (
            '<span style="color:#555;">Adjust PAL / JobTime and click Submit.</span>'
        )
        _calc_set_ops_status('', '#555')
        calc_update_view()

    def _calc_submit_xyz_browser_workflow(_button=None):
        mode = state.get('xyz_workflow_mode') or ''
        selected_path = state.get('xyz_workflow_source_path')
        if mode not in ('hyperpol_xtb', 'tadf_xtb') or selected_path is None:
            calc_xyz_workflow_status.value = (
                '<span style="color:#d32f2f;">No workflow prepared.</span>'
            )
            return

        selected_path = Path(selected_path)
        if not selected_path.exists():
            calc_xyz_workflow_status.value = (
                '<span style="color:#d32f2f;">Selected .xyz file no longer exists.</span>'
            )
            return
        if selected_path.suffix.lower() != '.xyz' or not _calc_is_single_structure_xyz(selected_path):
            calc_xyz_workflow_status.value = (
                '<span style="color:#d32f2f;">Select a single-structure .xyz file first.</span>'
            )
            return

        base_token = _calc_safe_job_token(selected_path.stem, 'xyz_job')
        if mode == 'hyperpol_xtb':
            folder_name = f'{base_token}_hyperpol_xtb'
            workflow_label = f'{base_token}_beta'
            submit_fn = ctx.backend.submit_hyperpol_xtb
        elif mode == 'tadf_xtb':
            folder_name = f'{base_token}_tadf_xtb'
            workflow_label = f'{base_token}_tadf'
            submit_fn = ctx.backend.submit_tadf_xtb
        else:
            calc_xyz_workflow_status.value = (
                f'<span style="color:#d32f2f;">Unknown workflow mode: {_html.escape(str(mode))}</span>'
            )
            return

        pal_used = max(1, int(calc_xyz_workflow_pal.value or CALC_BROWSER_WORKFLOW_PAL))
        time_limit = calc_xyz_workflow_time.value.strip() or CALC_BROWSER_WORKFLOW_TIMELIMIT
        target_dir = _calc_next_available_dir(selected_path.parent, folder_name)
        copied_xyz_name = f'{base_token}.xyz'
        copied_xyz_path = target_dir / copied_xyz_name
        try:
            target_dir.mkdir(parents=True, exist_ok=False)
            shutil.copy2(selected_path, copied_xyz_path)
        except Exception as exc:
            calc_xyz_workflow_status.value = (
                f'<span style="color:#d32f2f;">Failed to prepare '
                f'<code>{_html.escape(target_dir.name)}</code>: {_html.escape(str(exc))}</span>'
            )
            _calc_set_ops_status(
                f'Failed to prepare workflow folder <code>{_html.escape(target_dir.name)}</code>: '
                f'{_html.escape(str(exc))}',
                '#d32f2f',
            )
            return

        calc_xyz_workflow_status.value = (
            f'<span style="color:#1976d2;">Submitting <code>{_html.escape(mode)}</code> '
            f'from <code>{_html.escape(target_dir.name)}</code>...</span>'
        )
        _calc_set_ops_status(
            f'Submitting <code>{_html.escape(mode)}</code> from '
            f'<code>{_html.escape(target_dir.name)}</code>...',
            '#1976d2',
        )
        try:
            result = submit_fn(
                job_dir=target_dir,
                job_name=target_dir.name,
                xyz_file=copied_xyz_name,
                label=workflow_label,
                time_limit=time_limit,
                pal=pal_used,
                maxcore=CALC_BROWSER_WORKFLOW_MAXCORE,
                use_bfw=bool(calc_xyz_workflow_bfw.value),
            )
        except Exception as exc:
            calc_xyz_workflow_status.value = (
                f'<span style="color:#d32f2f;">Error submitting '
                f'<code>{_html.escape(mode)}</code>: {_html.escape(str(exc))}</span>'
            )
            _calc_set_ops_status(
                f'Error submitting <code>{_html.escape(mode)}</code>: {_html.escape(str(exc))}',
                '#d32f2f',
            )
            return

        calc_list_directory()
        if result.returncode == 0:
            submit_msg = (result.stdout or '').strip() or 'Submitted'
            calc_xyz_workflow_status.value = (
                f'<span style="color:#2e7d32;">Submitted <code>{_html.escape(mode)}</code> in '
                f'<code>{_html.escape(target_dir.name)}</code>: {_html.escape(submit_msg)}</span>'
            )
            _calc_set_ops_status(
                f'Submitted <code>{_html.escape(mode)}</code> in '
                f'<code>{_html.escape(target_dir.name)}</code>: {_html.escape(submit_msg)}',
                '#2e7d32',
            )
            calc_reset_xyz_workflow_state()
            _calc_reset_options_dropdown()
            calc_update_view()
        else:
            submit_msg = (result.stderr or result.stdout or 'Unknown error').strip()
            calc_xyz_workflow_status.value = (
                f'<span style="color:#d32f2f;">Failed to submit '
                f'<code>{_html.escape(mode)}</code> in <code>{_html.escape(target_dir.name)}</code>: '
                f'{_html.escape(submit_msg)}</span>'
            )
            _calc_set_ops_status(
                f'Failed to submit <code>{_html.escape(mode)}</code> in '
                f'<code>{_html.escape(target_dir.name)}</code>: {_html.escape(submit_msg)}',
                '#d32f2f',
            )

    def _calc_resolve_within_root(path):
        root = _calc_dir().resolve()
        resolved = path.resolve()
        try:
            resolved.relative_to(root)
        except ValueError as exc:
            raise ValueError('Path must remain inside explorer root.') from exc
        return resolved

    def _calc_extract_uploaded_zip(target_dir, filename, content):
        archive_name = str(filename or '').replace('\\', '/').split('/')[-1]
        archive_stem = Path(archive_name).stem or 'upload'
        extract_root = _calc_next_available_dir(target_dir, archive_stem)
        extract_root = _calc_resolve_within_root(extract_root)
        extract_root.mkdir(parents=False, exist_ok=False)

        extracted_files = []
        payload = bytes(content)
        with zipfile.ZipFile(io.BytesIO(payload)) as archive:
            file_entries = [info for info in archive.infolist() if not info.is_dir()]
            part_lists = []
            for info in file_entries:
                raw_name = str(info.filename or '').replace('\\', '/')
                if raw_name.startswith('/') or re.match(r'^[A-Za-z]:', raw_name):
                    raise ValueError(f'Unsafe ZIP entry: {info.filename}')
                parts = [part for part in raw_name.split('/') if part not in ('', '.')]
                if not parts or any(part == '..' for part in parts):
                    raise ValueError(f'Unsafe ZIP entry: {info.filename}')
                part_lists.append(parts)

            common_root = None
            if part_lists:
                first_root = part_lists[0][0]
                if all(parts and parts[0] == first_root for parts in part_lists):
                    common_root = first_root
            strip_common_root = bool(common_root) and any(len(parts) > 1 for parts in part_lists)

            for info, parts in zip(file_entries, part_lists):
                rel_parts = parts[1:] if strip_common_root and len(parts) > 1 else parts
                if not rel_parts:
                    rel_parts = [parts[-1]]
                destination = _calc_resolve_within_root(extract_root.joinpath(*rel_parts))
                destination.parent.mkdir(parents=True, exist_ok=True)
                destination = _calc_next_available_path(destination.parent, destination.name)
                with archive.open(info, 'r') as src, destination.open('wb') as dst:
                    shutil.copyfileobj(src, dst)
                extracted_files.append(destination)

        return extract_root, extracted_files

    def _calc_save_uploaded_entries(upload_entries, target_dir):
        saved_files = []
        extracted_archives = []
        for entry in upload_entries or ():
            name = entry['name'] if isinstance(entry, dict) else getattr(entry, 'name', '')
            content = entry['content'] if isinstance(entry, dict) else getattr(entry, 'content', b'')
            if not name:
                continue
            if str(name).lower().endswith('.zip'):
                extract_root, extracted_files = _calc_extract_uploaded_zip(target_dir, name, content)
                extracted_archives.append((extract_root, extracted_files))
                continue
            destination = _calc_next_available_path(target_dir, name)
            destination = _calc_resolve_within_root(destination)
            destination.parent.mkdir(parents=True, exist_ok=True)
            destination.write_bytes(bytes(content))
            saved_files.append(destination)
        return saved_files, extracted_archives

    def _calc_display_path(path):
        try:
            rel = Path(path).resolve().relative_to(_calc_dir().resolve())
        except Exception:
            return _html.escape(Path(path).name)
        rel_text = rel.as_posix()
        return '/' if rel_text in {'', '.'} else f'/{_html.escape(rel_text)}'

    def calc_on_hidden_upload(change):
        upload_entries = tuple(change.get('new') or ())
        if not upload_entries:
            return
        raw_target = str(calc_upload_target_input.value or '').strip()
        try:
            target_dir = _calc_current_dir() if not raw_target else _calc_parse_move_target_dir(raw_target)
            saved_files, extracted_archives = _calc_save_uploaded_entries(upload_entries, target_dir)
            summary_parts = []
            if saved_files:
                summary_parts.append(f'uploaded {len(saved_files)} file(s)')
            if extracted_archives:
                summary_parts.append(f'unpacked {len(extracted_archives)} ZIP archive(s)')
            if not summary_parts:
                summary_parts.append('no files received')
            _calc_set_ops_status(
                f'Explorer upload complete in <code>{_calc_display_path(target_dir)}</code>: '
                f'{"; ".join(summary_parts)}.',
                color='#2e7d32',
            )
            _calc_refresh_related_explorers()
        except Exception as exc:
            _calc_set_ops_status(
                f'Explorer upload failed: {_html.escape(str(exc))}',
                color='#d32f2f',
            )
        finally:
            calc_upload_target_input.value = ''
            try:
                calc_hidden_upload.value = ()
            except Exception:
                pass

    def _calc_normalize_upload_parts(raw_path):
        raw = str(raw_path or '').replace('\\', '/').strip('/')
        parts = []
        for part in raw.split('/'):
            piece = str(part or '').strip()
            if not piece or piece == '.':
                continue
            if piece == '..':
                raise ValueError('Upload path cannot contain "..".')
            parts.append(_calc_clean_item_name(piece))
        if not parts:
            raise ValueError('Upload path is empty.')
        return parts

    def _calc_store_uploaded_payload(target_dir, rel_parts, payload):
        parent_dir = target_dir
        if len(rel_parts) > 1:
            parent_dir = _calc_resolve_within_root(target_dir.joinpath(*rel_parts[:-1]))
            parent_dir.mkdir(parents=True, exist_ok=True)
        filename = rel_parts[-1]
        if filename.lower().endswith('.zip'):
            extract_root, extracted_files = _calc_extract_uploaded_zip(parent_dir, filename, payload)
            return {'kind': 'zip', 'path': extract_root, 'count': len(extracted_files)}
        destination = _calc_next_available_path(parent_dir, filename)
        destination = _calc_resolve_within_root(destination)
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_bytes(payload)
        return {'kind': 'file', 'path': destination, 'count': 1}

    def _calc_parse_staged_upload_target_dir(raw_target_rel):
        root_dir = _calc_dir().resolve()
        raw = str(raw_target_rel or '').replace('\\', '/').strip('/')
        if not raw:
            return root_dir

        rel_parts = []
        for part in raw.split('/'):
            piece = str(part or '').strip()
            if not piece or piece == '.':
                continue
            if piece == '..':
                raise ValueError('Upload target cannot contain "..".')
            rel_parts.append(_calc_clean_item_name(piece))

        if not rel_parts:
            return root_dir
        return _calc_resolve_within_root(root_dir.joinpath(*rel_parts))

    def _calc_load_staged_upload_manifest(batch_dir):
        manifest_path = Path(batch_dir) / 'manifest.json'
        if not manifest_path.is_file():
            return None
        return json.loads(manifest_path.read_text(encoding='utf-8'))

    def _calc_read_staged_upload_payload(batch_dir, staged_path):
        batch_path = Path(batch_dir).resolve()
        raw = str(staged_path or '').replace('\\', '/').strip('/')
        if not raw:
            raise ValueError('Staged upload path is empty.')
        payload_path = (batch_path / raw).resolve()
        try:
            payload_path.relative_to(batch_path)
        except Exception as exc:
            raise ValueError('Staged upload path escapes batch directory.') from exc
        if not payload_path.is_file():
            raise FileNotFoundError(f'Staged upload file not found: {payload_path.name}')
        return payload_path.read_bytes()

    def _calc_consume_staged_upload_batch(batch_dir):
        manifest = _calc_load_staged_upload_manifest(batch_dir)
        if not manifest:
            return None

        target_dir = _calc_parse_staged_upload_target_dir(manifest.get('target_dir_rel'))
        file_entries = manifest.get('files') or []
        if not file_entries:
            raise ValueError('Upload manifest does not contain any files.')

        saved_files = []
        unzipped_archives = []
        for item in file_entries:
            rel_parts = _calc_normalize_upload_parts(item.get('relative_path'))
            staged_path = item.get('staged_path') or (
                Path('payload').joinpath(*rel_parts).as_posix()
            )
            payload = _calc_read_staged_upload_payload(batch_dir, staged_path)
            saved = _calc_store_uploaded_payload(target_dir, rel_parts, payload)
            if saved['kind'] == 'zip':
                unzipped_archives.append(saved['path'])
            else:
                saved_files.append(saved['path'])

        return {
            'target_dir': target_dir,
            'saved_files': saved_files,
            'unzipped_archives': unzipped_archives,
        }

    def _calc_process_staged_uploads():
        staging_root = CALC_BROWSER_UPLOAD_STAGING_SCOPE_DIR
        if not staging_root.exists():
            return

        processed_batches = 0
        processed_files = 0
        processed_archives = 0
        last_target_dir = None
        errors = []

        try:
            batch_dirs = sorted(
                [
                    entry for entry in staging_root.iterdir()
                    if entry.is_dir() and not entry.name.endswith('.failed')
                ],
                key=lambda entry: entry.name,
            )
        except Exception:
            return

        for batch_dir in batch_dirs:
            if not (batch_dir / 'manifest.json').is_file():
                continue
            try:
                consumed = _calc_consume_staged_upload_batch(batch_dir)
                if not consumed:
                    continue
                processed_batches += 1
                processed_files += len(consumed['saved_files'])
                processed_archives += len(consumed['unzipped_archives'])
                last_target_dir = consumed['target_dir']
                shutil.rmtree(batch_dir, ignore_errors=True)
            except Exception as exc:
                errors.append(f'{batch_dir.name}: {exc}')
                try:
                    failed_dir = batch_dir.with_name(f'{batch_dir.name}.failed')
                    if not failed_dir.exists():
                        batch_dir.rename(failed_dir)
                except Exception:
                    pass

        if not processed_batches and not errors:
            return

        summary_parts = []
        if processed_files:
            summary_parts.append(f'uploaded {processed_files} file(s)')
        if processed_archives:
            summary_parts.append(f'unpacked {processed_archives} ZIP archive(s)')
        if errors:
            summary_parts.append(f'{len(errors)} batch error(s)')
        if not summary_parts:
            summary_parts.append('no files received')

        target_html = (
            f' in <code>{_calc_display_path(last_target_dir)}</code>'
            if last_target_dir is not None else ''
        )
        _calc_set_ops_status(
            f'Explorer upload complete{target_html}: {"; ".join(summary_parts)}.',
            color='#2e7d32' if not errors else '#d32f2f',
        )
        if processed_batches:
            _calc_refresh_related_explorers(include_self=False)

    def _calc_finalize_upload_batch(batch_id, target_dir):
        batch = upload_bridge_state['batches'].get(batch_id)
        if not batch:
            return
        summary_parts = []
        if batch['saved_files']:
            summary_parts.append(f'uploaded {len(batch["saved_files"])} file(s)')
        if batch['unzipped_archives']:
            summary_parts.append(f'unpacked {len(batch["unzipped_archives"])} ZIP archive(s)')
        if batch['errors']:
            summary_parts.append(f'{len(batch["errors"])} error(s)')
        if not summary_parts:
            summary_parts.append('no files received')
        status_color = '#2e7d32' if not batch['errors'] else '#d32f2f'
        _calc_set_ops_status(
            f'Explorer upload complete in <code>{_calc_display_path(target_dir)}</code>: '
            f'{"; ".join(summary_parts)}.',
            color=status_color,
        )
        _calc_refresh_related_explorers()
        upload_bridge_state['batches'].pop(batch_id, None)

    def calc_on_upload_bridge_seq(change):
        try:
            seq = int(change.get('new') or 0)
        except Exception:
            return
        if seq <= 0 or seq <= upload_bridge_state['last_seq']:
            return
        upload_bridge_state['last_seq'] = seq
        batch_id = ''
        try:
            meta = json.loads(calc_upload_meta_input.value or '{}')
            chunk_data = str(calc_upload_chunk_input.value or '')
            batch_id = str(meta.get('batch_id') or '').strip() or f'batch_{seq}'
            upload_id = str(meta.get('upload_id') or '').strip() or f'{batch_id}:file'
            batch_total = max(1, int(meta.get('batch_total') or 1))
            chunk_index = int(meta.get('chunk_index') or 0)
            chunk_total = max(1, int(meta.get('chunk_total') or 1))
            raw_target = str(meta.get('target') or '').strip()
            rel_parts = _calc_normalize_upload_parts(meta.get('relative_path') or meta.get('name'))
            target_dir = _calc_current_dir() if not raw_target else _calc_parse_move_target_dir(raw_target)

            batch = upload_bridge_state['batches'].setdefault(
                batch_id,
                {
                    'expected_files': batch_total,
                    'completed': set(),
                    'saved_files': [],
                    'unzipped_archives': [],
                    'errors': [],
                    'target_dir': target_dir,
                },
            )
            batch['expected_files'] = max(batch['expected_files'], batch_total)
            batch['target_dir'] = target_dir

            file_state = upload_bridge_state['files'].setdefault(
                upload_id,
                {
                    'batch_id': batch_id,
                    'target': raw_target,
                    'rel_parts': rel_parts,
                    'chunk_total': chunk_total,
                    'chunks': [None] * chunk_total,
                },
            )
            if file_state['chunk_total'] != chunk_total:
                raise ValueError('Upload chunk count changed during transfer.')
            if chunk_index < 0 or chunk_index >= chunk_total:
                raise ValueError('Invalid upload chunk index.')
            file_state['chunks'][chunk_index] = chunk_data

            if all(part is not None for part in file_state['chunks']):
                payload = base64.b64decode(''.join(file_state['chunks']))
                saved = _calc_store_uploaded_payload(target_dir, file_state['rel_parts'], payload)
                if saved['kind'] == 'zip':
                    batch['unzipped_archives'].append(saved['path'])
                else:
                    batch['saved_files'].append(saved['path'])
                batch['completed'].add(upload_id)
                upload_bridge_state['files'].pop(upload_id, None)
                if len(batch['completed']) >= batch['expected_files']:
                    _calc_finalize_upload_batch(batch_id, target_dir)
        except Exception as exc:
            if batch_id:
                batch = upload_bridge_state['batches'].setdefault(
                    batch_id,
                    {
                        'expected_files': 1,
                        'completed': set(),
                        'saved_files': [],
                        'unzipped_archives': [],
                        'errors': [],
                        'target_dir': _calc_current_dir(),
                    },
                )
                batch['errors'].append(str(exc))
            _calc_set_ops_status(
                f'Explorer upload failed: {_html.escape(str(exc))}',
                color='#d32f2f',
            )
        finally:
            calc_upload_ack_input.value = seq
            calc_upload_chunk_input.value = ''

    def calc_on_upload_trigger_click(_btn):
        """Handle upload chunk via button click (works in Voilà)."""
        upload_bridge_state['last_seq'] += 1
        seq = upload_bridge_state['last_seq']
        batch_id = ''
        try:
            meta = json.loads(calc_upload_meta_input.value or '{}')
            chunk_data = str(calc_upload_chunk_input.value or '')
            batch_id = str(meta.get('batch_id') or '').strip() or f'batch_{seq}'
            upload_id = str(meta.get('upload_id') or '').strip() or f'{batch_id}:file'
            batch_total = max(1, int(meta.get('batch_total') or 1))
            chunk_index = int(meta.get('chunk_index') or 0)
            chunk_total = max(1, int(meta.get('chunk_total') or 1))
            raw_target = str(meta.get('target') or '').strip()
            rel_parts = _calc_normalize_upload_parts(meta.get('relative_path') or meta.get('name'))
            target_dir = _calc_current_dir() if not raw_target else _calc_parse_move_target_dir(raw_target)

            batch = upload_bridge_state['batches'].setdefault(
                batch_id,
                {
                    'expected_files': batch_total,
                    'completed': set(),
                    'saved_files': [],
                    'unzipped_archives': [],
                    'errors': [],
                    'target_dir': target_dir,
                },
            )
            batch['expected_files'] = max(batch['expected_files'], batch_total)
            batch['target_dir'] = target_dir

            file_state = upload_bridge_state['files'].setdefault(
                upload_id,
                {
                    'batch_id': batch_id,
                    'target': raw_target,
                    'rel_parts': rel_parts,
                    'chunk_total': chunk_total,
                    'chunks': [None] * chunk_total,
                },
            )
            if file_state['chunk_total'] != chunk_total:
                raise ValueError('Upload chunk count changed during transfer.')
            if chunk_index < 0 or chunk_index >= chunk_total:
                raise ValueError('Invalid upload chunk index.')
            file_state['chunks'][chunk_index] = chunk_data

            if all(part is not None for part in file_state['chunks']):
                payload = base64.b64decode(''.join(file_state['chunks']))
                saved = _calc_store_uploaded_payload(target_dir, file_state['rel_parts'], payload)
                if saved['kind'] == 'zip':
                    batch['unzipped_archives'].append(saved['path'])
                else:
                    batch['saved_files'].append(saved['path'])
                batch['completed'].add(upload_id)
                upload_bridge_state['files'].pop(upload_id, None)
                if len(batch['completed']) >= batch['expected_files']:
                    _calc_finalize_upload_batch(batch_id, target_dir)
        except Exception as exc:
            if batch_id:
                batch = upload_bridge_state['batches'].setdefault(
                    batch_id,
                    {
                        'expected_files': 1,
                        'completed': set(),
                        'saved_files': [],
                        'unzipped_archives': [],
                        'errors': [],
                        'target_dir': _calc_current_dir(),
                    },
                )
                batch['errors'].append(str(exc))
            _calc_set_ops_status(
                f'Explorer upload failed: {_html.escape(str(exc))}',
                color='#d32f2f',
            )
        finally:
            calc_upload_ack_input.value = seq
            calc_upload_ack_label.value = str(seq)
            calc_upload_chunk_input.value = ''

    def _calc_is_relative_to(path, parent):
        try:
            path.relative_to(parent)
            return True
        except ValueError:
            return False

    def _calc_clean_item_name(raw_name):
        name = str(raw_name or '').strip()
        if not name:
            raise ValueError('Name cannot be empty.')
        if name in {'.', '..'}:
            raise ValueError('Name cannot be "." or "..".')
        if '/' in name or '\\' in name:
            raise ValueError('Name cannot contain path separators.')
        return name

    def _calc_next_available_folder_name(base_name='New Folder'):
        current_dir = _calc_current_dir()
        candidate = str(base_name)
        if not (current_dir / candidate).exists():
            return candidate
        idx = 2
        while True:
            candidate = f'{base_name} ({idx})'
            if not (current_dir / candidate).exists():
                return candidate
            idx += 1

    def _calc_selected_labels():
        return [label for label in calc_file_list.value if label and not label.startswith('(')]

    def _calc_take_forced_action_source():
        forced = str(calc_action_source_input.value or '').strip()
        calc_action_source_input.value = ''
        if not forced:
            return None
        try:
            name = _calc_clean_item_name(forced)
        except Exception:
            return None
        src = _calc_current_dir() / name
        return src if src.exists() else None

    def _calc_collect_move_sources():
        forced_src = _calc_take_forced_action_source()
        if forced_src is not None:
            return [forced_src]
        labels = _calc_selected_labels()
        current_dir = _calc_current_dir()
        sources = []
        for label in labels:
            name = _calc_label_to_name(label)
            if not name:
                continue
            src = current_dir / name
            if src.exists():
                sources.append(src)
        if sources:
            return sources
        # Fallback: allow moving current folder when nothing is selected.
        if state['current_path']:
            cur = _calc_current_dir()
            if cur.exists() and cur.resolve() != _calc_dir().resolve():
                return [cur]
        return []

    def _calc_collect_selected_sources_only():
        labels = _calc_selected_labels()
        current_dir = _calc_current_dir()
        sources = []
        for label in labels:
            name = _calc_label_to_name(label)
            if not name:
                continue
            src = current_dir / name
            if src.exists():
                sources.append(src)
        return sources

    def _calc_collect_rename_source():
        forced_src = _calc_take_forced_action_source()
        if forced_src is not None:
            return forced_src
        labels = _calc_selected_labels()
        if len(labels) > 1:
            raise ValueError('Select exactly one item for rename.')
        if labels:
            name = _calc_label_to_name(labels[0])
            if not name:
                raise ValueError('Invalid selection.')
            src = _calc_current_dir() / name
            if not src.exists():
                raise ValueError('Selected item not found.')
            return src
        # Fallback: rename current folder when nothing is selected.
        if state['current_path']:
            src = _calc_current_dir()
            if src.exists() and src.resolve() != _calc_dir().resolve():
                return src
        raise ValueError('Select one item or navigate into a folder to rename it.')

    def _calc_collect_duplicate_source():
        forced_src = _calc_take_forced_action_source()
        if forced_src is not None:
            src = forced_src
        else:
            labels = _calc_selected_labels()
            if len(labels) > 1:
                raise ValueError('Select exactly one folder to duplicate.')
            if labels:
                name = _calc_label_to_name(labels[0])
                if not name:
                    raise ValueError('Invalid selection.')
                src = _calc_current_dir() / name
                if not src.exists():
                    raise ValueError('Selected folder not found.')
            elif state['current_path']:
                src = _calc_current_dir()
                if not src.exists() or src.resolve() == _calc_dir().resolve():
                    raise ValueError('Select one folder or open a folder to duplicate it.')
            else:
                raise ValueError('Select one folder or open a folder to duplicate it.')
        src_resolved = _calc_resolve_within_root(src)
        if not src_resolved.is_dir():
            raise ValueError('Duplicate works only for folders.')
        return src_resolved

    def _calc_refresh_listing_preserve_filter():
        saved_filter = calc_folder_search.value
        calc_list_directory()
        if saved_filter:
            calc_folder_search.value = saved_filter
            calc_filter_file_list()

    def _calc_shared_explorer_state():
        shared = getattr(ctx, 'shared_explorer_state', None)
        if not isinstance(shared, dict):
            shared = {'refresh_hooks': {}, 'refresh_running': False}
            ctx.shared_explorer_state = shared
        hooks = shared.get('refresh_hooks')
        if not isinstance(hooks, dict):
            hooks = {}
            shared['refresh_hooks'] = hooks
        if 'refresh_running' not in shared:
            shared['refresh_running'] = False
        return shared

    def _calc_register_explorer_refresh():
        shared = _calc_shared_explorer_state()
        shared['refresh_hooks'][calc_scope_id] = _calc_refresh_listing_preserve_filter

    def _calc_refresh_related_explorers(include_self=True):
        shared = _calc_shared_explorer_state()
        if shared.get('refresh_running'):
            return
        shared['refresh_running'] = True
        try:
            hooks = list((shared.get('refresh_hooks') or {}).items())
            for scope_key, callback in hooks:
                if not include_self and scope_key == calc_scope_id:
                    continue
                try:
                    callback()
                except Exception:
                    continue
        finally:
            shared['refresh_running'] = False

    def _calc_transfer_jobs_dir():
        return ensure_jobs_dir()

    def _calc_settings_path():
        return get_settings_path()

    def _calc_load_transfer_config():
        try:
            return load_transfer_settings()
        except Exception as exc:
            ctx.select_tab('Settings')
            _calc_set_ops_status(
                (
                    f'Settings invalid in <code>{_html.escape(str(_calc_settings_path()))}</code>: '
                    f'{_html.escape(str(exc))}'
                ),
                color='#d32f2f',
            )
            return False

    def _calc_transfer_status_color(status):
        mapping = {
            'queued': '#1976d2',
            'running': '#1976d2',
            'retrying': '#ef6c00',
            'success': '#2e7d32',
            'warning': '#ef6c00',
            'failed': '#d32f2f',
        }
        return mapping.get(str(status or '').lower(), '#555')

    def _calc_update_transfer_jobs_visibility():
        jobs_visible = calc_transfer_jobs_panel.layout.display != 'none'
        calc_transfer_jobs_btn.button_style = 'primary' if jobs_visible else ''
        if jobs_visible:
            calc_left.add_class('calc-transfer-jobs-mode')
        else:
            calc_left.remove_class('calc-transfer-jobs-mode')
        if jobs_visible:
            calc_rename_prompt_row.layout.display = 'none'
            calc_duplicate_prompt_row.layout.display = 'none'
            calc_filter_sort_row.layout.display = 'none'
            calc_file_list.layout.display = 'none'
        else:
            calc_filter_sort_row.layout.display = 'flex'
            calc_file_list.layout.display = ''

    def _calc_transfer_items_html(entry, limit=4):
        raw_sources = entry.get('sources', []) or []
        rendered = []
        for raw_source in raw_sources[:limit]:
            source_text = str(raw_source or '').strip()
            if not source_text:
                continue
            source_path = Path(source_text)
            label = source_path.name or source_text
            rendered.append(
                f'<code title="{_html.escape(source_text)}">{_html.escape(label)}</code>'
            )
        if not rendered:
            return 'n/a'
        remaining = len(raw_sources) - len(rendered)
        summary = ', '.join(rendered)
        if remaining > 0:
            summary = f'{summary} + {remaining} more'
        return summary

    def _calc_render_transfer_jobs(limit=8):
        try:
            entries = list_transfer_jobs(jobs_dir=_calc_transfer_jobs_dir(), limit=limit)
        except Exception as exc:
            calc_transfer_jobs_html.value = (
                f'<span style="color:#d32f2f;">Could not load transfer jobs: '
                f'{_html.escape(str(exc))}</span>'
            )
            return

        if not entries:
            calc_transfer_jobs_html.value = (
                '<span style="color:#555;">No background transfer jobs yet.</span>'
            )
            return

        blocks = []
        for entry in entries:
            status = str(entry.get('status') or 'unknown')
            color = _calc_transfer_status_color(status)
            direction = str(entry.get('direction') or 'push').lower()
            remote = (
                f'{entry.get("user", "?")}@{entry.get("host", "?")}:'
                f'{entry.get("remote_path", "?")}'
            )
            local_target = str(entry.get('local_target') or '').strip()
            endpoint = remote
            if direction == 'pull' and local_target:
                endpoint = f'{remote} -> {local_target}'
            summary = str(entry.get('last_summary') or entry.get('last_error') or '').strip()
            retry_note = ''
            retry_in = entry.get('retry_in_seconds', 0) or 0
            if status == 'retrying' and retry_in:
                retry_note = f' Retry in about {int(retry_in)} s.'
            attempts = entry.get('attempt', 0)
            max_retries = entry.get('max_retries', 0)
            log_path = entry.get('log_path', '')
            items_html = _calc_transfer_items_html(entry)
            updated = str(entry.get('updated_at') or '').replace('T', ' ')
            if updated.endswith('+00:00'):
                updated = updated[:-6] + ' UTC'
            blocks.append(
                '<div style="border:1px solid #d9dee3;border-radius:6px;padding:8px 10px;'
                'margin:0 0 8px 0;background:#fafbfc;">'
                f'<div><b>{_html.escape(entry.get("job_id", "job"))}</b> '
                f'<span style="color:{color};font-weight:600;">{_html.escape(status.upper())}</span></div>'
                f'<div><code>{_html.escape(endpoint)}</code></div>'
                f'<div>{int(entry.get("source_count", 0) or 0)} item(s), '
                f'attempt {int(attempts or 0)}/{int(max_retries or 0) + 1}</div>'
                f'<div>Items: {items_html}</div>'
                f'<div>{_html.escape(summary or "No status message.")}{_html.escape(retry_note)}</div>'
                f'<div style="color:#555;">Updated: {_html.escape(updated or "-")}</div>'
                f'<div style="color:#555;">Log: <code>{_html.escape(str(log_path))}</code></div>'
                '</div>'
            )
        calc_transfer_jobs_html.value = ''.join(blocks)

    def _calc_set_ssh_transfer_running(is_running):
        state['ssh_transfer_running'] = bool(is_running)
        _calc_update_explorer_action_state()

    def _calc_hide_rename_prompt():
        state['rename_source_path'] = ''
        calc_rename_prompt_input.value = ''
        calc_rename_prompt_row.layout.display = 'none'
        state['duplicate_source_path'] = ''
        calc_duplicate_prompt_input.value = ''
        calc_duplicate_prompt_row.layout.display = 'none'
        calc_new_folder_prompt_input.value = ''
        calc_new_folder_prompt_row.layout.display = 'none'

    def _calc_clipboard_store():
        store = getattr(ctx, 'shared_clipboard', None)
        if not isinstance(store, dict):
            store = {
                'paths': list(getattr(ctx, 'clipboard_paths', []) or []),
                'mode': str(getattr(ctx, 'clipboard_mode', '') or ''),
            }
            ctx.shared_clipboard = store
        if 'paths' not in store:
            store['paths'] = []
        if 'mode' not in store:
            store['mode'] = ''
        return store

    def _calc_clipboard_paths():
        store = _calc_clipboard_store()
        raw_paths = store.get('paths', []) or []
        valid = []
        for raw in raw_paths:
            try:
                p = Path(raw).resolve()
            except Exception:
                continue
            if p.exists():
                valid.append(p)
        normalized = [str(p) for p in valid]
        store['paths'] = normalized
        ctx.clipboard_paths = normalized
        return valid

    def _calc_clipboard_mode():
        store = _calc_clipboard_store()
        mode = str(store.get('mode', '') or '')
        ctx.clipboard_mode = mode
        return mode

    def _calc_update_explorer_action_state():
        _calc_clipboard_paths()
        if _remote_archive_enabled:
            has_selection = bool(_calc_collect_selected_sources_only())
            calc_ssh_transfer_btn.disabled = (
                state.get('ssh_transfer_running', False) or not has_selection
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
                        {VIEWER_MOUSE_PATCH_JS}
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
            calc_xyz_workflow_toolbar.layout.display = 'none'
            calc_nmr_panel.layout.display = 'none'
            calc_censo_nmr_panel.layout.display = 'none'
        else:
            calc_preselect_container.layout.display = 'none'
            # Restore toolbars based on current mode.
            if state['recalc_active']:
                calc_content_toolbar.layout.display = 'none'
                calc_recalc_toolbar.layout.display = 'flex'
                calc_xyz_workflow_toolbar.layout.display = 'none'
            elif state['xyz_workflow_active']:
                calc_content_toolbar.layout.display = 'none'
                calc_recalc_toolbar.layout.display = 'none'
                calc_xyz_workflow_toolbar.layout.display = 'flex'
            elif calc_view_toggle.value:
                calc_content_toolbar.layout.display = 'none'
                calc_recalc_toolbar.layout.display = 'none'
                calc_xyz_workflow_toolbar.layout.display = 'none'
            else:
                calc_content_toolbar.layout.display = 'flex'
                calc_recalc_toolbar.layout.display = 'none'
                calc_xyz_workflow_toolbar.layout.display = 'none'
            calc_nmr_panel.layout.display = 'none'
            calc_censo_nmr_panel.layout.display = 'none'

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

    def _calc_build_png_filename():
        selected_path = _calc_get_selected_path()
        if selected_path:
            base = selected_path.stem or 'viewer'
        else:
            base = 'viewer'
        base = re.sub(r'[^A-Za-z0-9._-]+', '_', str(base)).strip('._') or 'viewer'
        if len(state.get('xyz_frames') or []) > 1:
            frame_no = int(state.get('xyz_current_frame', [0])[0]) + 1
            return f'{base}_frame_{frame_no:04d}.png'
        return f'{base}.png'

    def _calc_set_png_button_mode(main=False):
        calc_view_png_btn.layout.display = 'inline-flex' if main else 'none'
        calc_xyz_png_row.layout.display = 'none'
        calc_coord_png_row.layout.display = 'none'

    def _calc_trigger_png_download():
        scope_key_json = json.dumps(calc_scope_id)
        filename_json = json.dumps(_calc_build_png_filename())
        _run_js(
            f"""
            (function() {{
                var scopeKey = {scope_key_json};
                var filename = {filename_json};
                var viewerMap = window._calcMolViewerByScope || {{}};
                var trajMap = window._calcTrajViewerByScope || {{}};
                var viewer = trajMap[scopeKey] || viewerMap[scopeKey] || null;
                if (!viewer || !window.__delfinDownloadViewerPng) return;
                window.__delfinDownloadViewerPng(viewer, {{
                    filename: filename,
                    scale: {CALC_VIEWER_PNG_SCALE}
                }});
            }})();
            """
        )

    def _calc_traj_can_play():
        n_frames = len(state.get('xyz_frames') or [])
        return 1 < n_frames <= CALC_XYZ_LARGE_TRAJ_FRAMES

    def _calc_update_loop_button_style():
        calc_xyz_loop_checkbox.button_style = 'info' if calc_xyz_loop_checkbox.value else ''

    def _calc_set_play_button_state(active, sync_value=True):
        active = bool(active)
        state['traj_playing'] = active
        calc_xyz_play_btn.description = 'Stop' if active else 'Play'
        calc_xyz_play_btn.icon = 'stop' if active else 'play'
        calc_xyz_play_btn.button_style = 'danger' if active else 'success'
        if sync_value and calc_xyz_play_btn.value != active:
            state['traj_play_toggle_guard'] = True
            try:
                calc_xyz_play_btn.value = active
            finally:
                state['traj_play_toggle_guard'] = False

    def _calc_stop_xyz_playback(update_button=True):
        stop_event = state.get('traj_play_stop_event')
        state['traj_play_stop_event'] = None
        if isinstance(stop_event, threading.Event):
            stop_event.set()
        scope_key_json = json.dumps(calc_scope_id)
        _run_js(
            f"""
            (function() {{
                var scopeKey = {scope_key_json};
                if (window._calcTrajPlayTimerByScope && window._calcTrajPlayTimerByScope[scopeKey]) {{
                    clearInterval(window._calcTrajPlayTimerByScope[scopeKey]);
                    delete window._calcTrajPlayTimerByScope[scopeKey];
                }}
            }})();
            """
        )
        if update_button:
            _calc_set_play_button_state(False, sync_value=True)

    def _calc_update_xyz_traj_control_state():
        n_frames = len(state.get('xyz_frames') or [])
        can_play = _calc_traj_can_play()
        xyz_controls_visible = str(calc_xyz_controls.layout.display or 'none') != 'none'
        coord_controls_visible = str(calc_coord_controls.layout.display or 'none') != 'none'
        show_tray = xyz_controls_visible or coord_controls_visible
        calc_xyz_loop_checkbox.disabled = not can_play
        calc_xyz_fps_input.disabled = not can_play
        calc_xyz_play_btn.disabled = not can_play
        calc_xyz_tray_controls.layout.display = 'flex' if show_tray else 'none'
        _calc_update_loop_button_style()
        if not can_play:
            _calc_stop_xyz_playback(update_button=True)

    def _calc_apply_traj_style():
        style_js = DEFAULT_3DMOL_STYLE_JS
        scope_key_json = json.dumps(calc_scope_id)
        _run_js(
            f"""
            (function() {{
                var scopeKey = {scope_key_json};
                var viewer = null;
                if (window._calcTrajViewerByScope && window._calcTrajViewerByScope[scopeKey]) {{
                    viewer = window._calcTrajViewerByScope[scopeKey];
                }} else if (window._calcMolViewerByScope && window._calcMolViewerByScope[scopeKey]) {{
                    viewer = window._calcMolViewerByScope[scopeKey];
                }}
                if (!viewer) return;
                try {{
                    viewer.setStyle({{}}, {style_js});
                    viewer.render();
                }} catch (_e) {{}}
            }})();
            """
        )

    def _calc_render_traj_frame(frame_idx):
        scope_key_json = json.dumps(calc_scope_id)
        frame_idx = int(frame_idx)
        _run_js(
            f"""
            (function() {{
                var scopeKey = {scope_key_json};
                var tries = 0;
                function applyFrame() {{
                    var viewer = null;
                    if (window._calcTrajViewerByScope && window._calcTrajViewerByScope[scopeKey]) {{
                        viewer = window._calcTrajViewerByScope[scopeKey];
                    }} else if (window._calcMolViewerByScope && window._calcMolViewerByScope[scopeKey]) {{
                        viewer = window._calcMolViewerByScope[scopeKey];
                    }}
                    if (!viewer) {{
                        tries += 1;
                        if (tries < 24) setTimeout(applyFrame, 50);
                        return;
                    }}
                    try {{
                        viewer.setFrame({frame_idx});
                        viewer.render();
                    }} catch (_e) {{}}
                }}
                setTimeout(applyFrame, 0);
            }})();
            """
        )

    def _calc_start_xyz_playback():
        if not _calc_traj_can_play():
            _calc_set_play_button_state(False, sync_value=True)
            return

        if not state.get('traj_viewer_ready'):
            calc_update_xyz_viewer(initial_load=True)

        _calc_stop_xyz_playback(update_button=False)
        _calc_set_play_button_state(True, sync_value=False)
        frames = state.get('xyz_frames') or []
        frame_count = len(frames)
        if frame_count <= 1:
            _calc_set_play_button_state(False, sync_value=True)
            return
        try:
            fps = int(calc_xyz_fps_input.value)
        except Exception:
            fps = CALC_XYZ_PLAY_FPS_DEFAULT
        fps = max(CALC_XYZ_PLAY_FPS_MIN, min(CALC_XYZ_PLAY_FPS_MAX, fps))
        delay_ms = max(16, int(round(1000.0 / float(fps))))
        loop_enabled = bool(calc_xyz_loop_checkbox.value)
        start_frame = int(state['xyz_current_frame'][0]) + 1
        scope_key_json = json.dumps(calc_scope_id)
        _run_js(
            f"""
            (function() {{
                var scopeKey = {scope_key_json};
                var frameCount = {int(frame_count)};
                var loopEnabled = {str(loop_enabled).lower()};
                var delayMs = {int(delay_ms)};
                var startFrame = {int(start_frame)};
                window._calcTrajPlayTimerByScope = window._calcTrajPlayTimerByScope || {{}};
                if (window._calcTrajPlayTimerByScope[scopeKey]) {{
                    clearInterval(window._calcTrajPlayTimerByScope[scopeKey]);
                    delete window._calcTrajPlayTimerByScope[scopeKey];
                }}
                function getViewer() {{
                    if (window._calcTrajViewerByScope && window._calcTrajViewerByScope[scopeKey]) {{
                        return window._calcTrajViewerByScope[scopeKey];
                    }}
                    if (window._calcMolViewerByScope && window._calcMolViewerByScope[scopeKey]) {{
                        return window._calcMolViewerByScope[scopeKey];
                    }}
                    return null;
                }}
                var scopeRoot = document.querySelector('.{calc_scope_id}');
                if (!scopeRoot) return;
                var frameWidget = scopeRoot.querySelector('.calc-xyz-frame-input');
                var frameInput = frameWidget ? frameWidget.querySelector('input') : null;
                var playButton = scopeRoot.querySelector('.calc-xyz-play-btn button');
                if (!frameInput) return;
                var frameModel = null;
                (function initFrameModel() {{
                    try {{
                        var modelId = frameWidget
                            ? (frameWidget.getAttribute('data-widget-id')
                                || frameWidget.getAttribute('data-model-id'))
                            : null;
                        if (!modelId) return;
                        var mgr = null;
                        if (window.Jupyter
                            && window.Jupyter.notebook
                            && window.Jupyter.notebook.kernel
                            && window.Jupyter.notebook.kernel.widget_manager) {{
                            mgr = window.Jupyter.notebook.kernel.widget_manager;
                        }} else if (window.jupyterWidgetManager) {{
                            mgr = window.jupyterWidgetManager;
                        }}
                        if (!mgr || typeof mgr.get_model !== 'function') return;
                        var modelPromise = mgr.get_model(modelId);
                        if (modelPromise && typeof modelPromise.then === 'function') {{
                            modelPromise.then(function(model) {{
                                frameModel = model || null;
                            }});
                        }}
                    }} catch (_e) {{}}
                }})();

                function syncFrameValue(next) {{
                    frameInput.value = String(next);
                    frameInput.dispatchEvent(new Event('input', {{bubbles: true}}));
                    frameInput.dispatchEvent(new Event('change', {{bubbles: true}}));
                    if (frameModel) {{
                        try {{
                            frameModel.set('value', next);
                            frameModel.save_changes();
                        }} catch (_e) {{}}
                    }}
                }}

                var current = parseInt(frameInput.value || String(startFrame), 10);
                if (!isFinite(current) || current < 1 || current > frameCount) {{
                    current = Math.max(1, Math.min(frameCount, startFrame));
                    syncFrameValue(current);
                }}
                var initViewer = getViewer();
                if (initViewer) {{
                    try {{
                        initViewer.setFrame(current - 1);
                        initViewer.render();
                    }} catch (_e) {{}}
                }}

                var timer = setInterval(function() {{
                    var nowVal = parseInt(frameInput.value || String(current), 10);
                    if (!isFinite(nowVal) || nowVal < 1 || nowVal > frameCount) {{
                        nowVal = current;
                    }}
                    var next = nowVal + 1;
                    if (next > frameCount) {{
                        if (loopEnabled) {{
                            next = 1;
                        }} else {{
                            clearInterval(timer);
                            delete window._calcTrajPlayTimerByScope[scopeKey];
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
                window._calcTrajPlayTimerByScope[scopeKey] = timer;
            }})();
            """
        )

    def _calc_copy_to_clipboard(text, label='content'):
        text_payload = json.dumps(text)
        label_payload = json.dumps(label)
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
            "}catch(_e){"
            "return false;"
            "}"
            "}"
            "function _done(method){"
            "console.log('Copied ' + label + ' via ' + method);"
            "}"
            "if(navigator.clipboard && navigator.clipboard.writeText){"
            "navigator.clipboard.writeText(text)"
            ".then(function(){_done('navigator.clipboard');})"
            ".catch(function(err){"
            "console.warn('Clipboard API failed for ' + label + ':', err);"
            "if(_legacyCopy()){_done('execCommand');return;}"
            "_manualPrompt();"
            "});"
            "}else{"
            "if(_legacyCopy()){_done('execCommand');return;}"
            "_manualPrompt();"
            "}"
            "})();"
        )

    def _calc_clear_main_viewer_state(reset_view_state=False):
        calc_mol_viewer.clear_output()
        _calc_set_png_button_mode(main=False)
        scope_key_json = json.dumps(calc_scope_id)
        clear_views_flag = 'true' if reset_view_state else 'false'
        current_view_scope_json = json.dumps(
            f"{calc_scope_id}:{state.get('current_path') or '/'}"
        )
        _run_js(
            f"""
            (function() {{
                var scopeKey = {scope_key_json};
                var currentViewScope = {current_view_scope_json};
                window._calcMolViewStateByScope = window._calcMolViewStateByScope || {{}};
                window._calcMolViewScopeKeyByScope = window._calcMolViewScopeKeyByScope || {{}};
                if (!{clear_views_flag}) {{
                    var currentViewer =
                        (window._calcMolViewerByScope && window._calcMolViewerByScope[scopeKey])
                        || (window._calcTrajViewerByScope && window._calcTrajViewerByScope[scopeKey])
                        || null;
                    var previousScope =
                        window._calcMolViewScopeKeyByScope[scopeKey] || currentViewScope;
                    if (currentViewer && typeof currentViewer.getView === 'function') {{
                        try {{
                            window._calcMolViewStateByScope[previousScope] = currentViewer.getView();
                        }} catch (_e) {{}}
                    }}
                }}
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
                /* Compute size from actual free space in the right panel */
                var rightPanel = scopeRoot ? scopeRoot.querySelector('.calc-right') : null;
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
                                if (!style || style.display === 'none' || style.visibility === 'hidden') {{
                                    continue;
                                }}
                                var childRect = child.getBoundingClientRect();
                                if (childRect.height > 0) reservedBelow += childRect.height;
                            }}
                        }}
                        var availH = rightRect.bottom - mvRect.top - reservedBelow - 10;
                        var row = mv.closest('.calc-mol-view-row');
                        var rowRect = row ? row.getBoundingClientRect() : rightRect;
                        var tray = scopeRoot ? scopeRoot.querySelector('.calc-xyz-tray-controls') : null;
                        var trayStyle = tray ? window.getComputedStyle(tray) : null;
                        var trayVisible = !!(tray && trayStyle && trayStyle.display !== 'none');
                        var trayWidth = trayVisible ? tray.getBoundingClientRect().width : 0;
                        var availW = Math.max(120, rowRect.width - trayWidth - 16);
                        var h = Math.floor(availH * {CALC_MOL_DYNAMIC_SCALE});
                        var w = Math.floor(Math.min(h * 1.2, availW));
                        if (h >= 80 && w >= 120) {{
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
                {VIEWER_MOUSE_PATCH_JS}
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
        _calc_set_png_button_mode(main=True)
        _calc_set_png_button_mode(main=True)

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
        labels = _calc_selected_labels()
        if not labels:
            return False
        name = _calc_label_to_name(labels[0])
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
        labels = _calc_selected_labels()
        if not labels:
            return None
        name = _calc_label_to_name(labels[0])
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

    def _calc_xyz_batch_default_filename():
        selected_path = _calc_selected_item_path()
        if selected_path and selected_path.suffix.lower() == '.xyz':
            base = selected_path.stem
        elif state.get('xyz_batch_path'):
            base = Path(state['xyz_batch_path']).name
        elif state.get('current_path'):
            base = Path(state['current_path']).name
        else:
            base = 'xyz_batch'
        safe = ''.join(c for c in str(base) if c.isalnum() or c in ('_', '-'))
        if not safe:
            safe = 'xyz_batch'
        if not safe.endswith('_batch'):
            safe = f'{safe}_batch'
        return f'{safe}.txt'

    def _calc_xyz_batch_clean_xyz(raw_xyz):
        text = (raw_xyz or '').strip()
        if not text:
            return ''
        try:
            frames = parse_xyz_frames(text)
        except Exception:
            frames = []
        if frames:
            _comment, xyz_block, _n_atoms = frames[0]
            return '\n'.join(line.rstrip() for line in xyz_block.splitlines() if line.strip()).strip()
        lines = text.splitlines()
        if len(lines) >= 3:
            try:
                int(lines[0].strip())
                return '\n'.join(line.rstrip() for line in lines[2:] if line.strip()).strip()
            except ValueError:
                pass
        return '\n'.join(line.rstrip() for line in lines if line.strip()).strip()

    def _calc_xyz_batch_resolve_source(item_path):
        if item_path.is_file():
            if item_path.suffix.lower() == '.xyz':
                return item_path, None
            return None, f'{item_path.name}: not an .xyz file'
        if item_path.is_dir():
            named_xyz = item_path / f'{item_path.name}.xyz'
            if named_xyz.exists() and named_xyz.is_file():
                return named_xyz, None
            try:
                xyz_candidates = sorted(
                    [p for p in item_path.iterdir() if p.is_file() and p.suffix.lower() == '.xyz'],
                    key=lambda p: p.name.lower(),
                )
            except Exception:
                return None, f'{item_path.name}: cannot scan folder'
            if xyz_candidates:
                return xyz_candidates[0], None
            return None, f'{item_path.name}: no .xyz file found'
        return None, f'{item_path.name}: unsupported selection'

    def _calc_xyz_batch_dir():
        rel_path = str(state.get('xyz_batch_path') or '').strip().strip('/')
        if not rel_path:
            return _calc_dir()
        try:
            candidate = _calc_resolve_within_root(_calc_dir() / rel_path)
        except Exception:
            state['xyz_batch_path'] = ''
            return _calc_dir()
        if candidate.exists() and candidate.is_dir():
            return candidate
        state['xyz_batch_path'] = ''
        return _calc_dir()

    def _calc_set_xyz_batch_dir(target_dir):
        try:
            resolved = _calc_resolve_within_root(Path(target_dir))
        except Exception:
            return False
        if not resolved.exists() or not resolved.is_dir():
            return False
        state['xyz_batch_path'] = _calc_state_rel_path(resolved)
        return True

    def _calc_collect_xyz_batch_labels(base_dir):
        labels = []
        marked = set(state.get('xyz_batch_marked_paths') or [])
        try:
            entries = sorted(
                base_dir.iterdir(),
                key=lambda x: (not x.is_dir(), x.name.lower()),
            )
        except Exception:
            return labels

        for entry in entries:
            rel_value = _calc_state_rel_path(entry)
            mark = '✔' if rel_value in marked else '□'
            if entry.is_dir():
                labels.append((f'{mark} 📂 {entry.name}', rel_value))
            elif entry.is_file() and entry.suffix.lower() == '.xyz':
                labels.append((f'{mark} 🔬 {entry.name}', rel_value))
        return labels

    def _calc_refresh_xyz_batch_selector():
        root_dir = _calc_xyz_batch_dir()
        candidates = _calc_collect_xyz_batch_labels(root_dir)
        previous_selection = set(calc_xyz_batch_select.value or ())
        # Save scroll position before re-render, restore after.
        _run_js(
            "var el=document.querySelector('.calc-xyz-batch-select select');"
            "if(el) window._batchScrollTop=el.scrollTop;"
        )
        # Force a real widget rerender so changed labels (⬜ -> ✅) become visible.
        calc_xyz_batch_select.options = ()
        calc_xyz_batch_select.value = ()
        calc_xyz_batch_select.options = tuple(candidates)
        kept = tuple(value for _label, value in candidates if value in previous_selection)
        selected_value = ''
        selected_path = _calc_selected_item_path()
        if selected_path is not None:
            try:
                if selected_path.parent.resolve() == root_dir.resolve():
                    selected_value = _calc_state_rel_path(selected_path)
            except Exception:
                selected_value = ''
        if kept:
            calc_xyz_batch_select.value = kept
        elif selected_value and any(value == selected_value for _label, value in candidates):
            calc_xyz_batch_select.value = (selected_value,)
        else:
            calc_xyz_batch_select.value = ()
        current_rel = state.get('xyz_batch_path') or ''
        display_rel = f'/{current_rel}' if current_rel else '/'
        marked_count = len(state.get('xyz_batch_marked_paths') or [])
        calc_xyz_batch_root_info.value = (
            f'<span style="color:#555;">Batch folder:</span> '
            f'<code>{_html.escape(display_rel)}</code> '
            f'<span style="color:#555;">| Haken:</span> <code>{marked_count}</code>'
        )
        calc_xyz_batch_up_btn.disabled = not bool(current_rel)
        calc_xyz_batch_root_btn.disabled = not bool(current_rel)
        calc_xyz_batch_select_all_btn.disabled = not bool(candidates)
        if not calc_xyz_batch_filename.value.strip():
            calc_xyz_batch_filename.value = _calc_xyz_batch_default_filename()
        # Restore scroll position after widget re-render
        _run_js(
            "setTimeout(function(){"
            "var el=document.querySelector('.calc-xyz-batch-select select');"
            "if(el && window._batchScrollTop!=null) el.scrollTop=window._batchScrollTop;"
            "}, 50);"
        )

    def _calc_show_xyz_batch_panel(show):
        if show:
            _calc_refresh_xyz_batch_selector()
            calc_xyz_batch_panel.layout.display = 'block'
        else:
            calc_xyz_batch_panel.layout.display = 'none'
        _run_js(
            f"setTimeout(function(){{"
            f" if(window['{calc_resize_mol_fn}']) window['{calc_resize_mol_fn}']();"
            f"}}, 80);"
        )

    def _calc_print_mode_parse_modes(raw_modes):
        tokens = [tok.strip() for tok in re.split(r'[\s,]+', str(raw_modes or '').strip()) if tok.strip()]
        if not tokens:
            return None, 'Please enter at least one mode (e.g. 6 7 8).'
        modes = []
        for token in tokens:
            try:
                mode_idx = int(token)
            except ValueError:
                return None, f'Invalid mode: {token}'
            if mode_idx <= 0:
                return None, 'Mode indices must be positive integers.'
            modes.append(str(mode_idx))
        return modes, ''

    def _calc_print_mode_target_paths():
        out_path = _calc_get_selected_path()
        if out_path is None or not out_path.exists():
            return None, None, 'Select a .out file first.'
        if out_path.suffix.lower() != '.out':
            return None, None, 'Selected file is not a .out file.'
        hess_path = out_path.with_suffix('.hess')
        if not hess_path.exists():
            return out_path, None, f'Missing Hessian file: {hess_path.name}'
        return out_path, hess_path, ''

    def _calc_refresh_print_mode_info():
        out_path = _calc_get_selected_path()
        if out_path and out_path.suffix.lower() == '.out':
            hess_path = out_path.with_suffix('.hess')
            if hess_path.exists():
                hess_html = f'<span style="color:#2e7d32;">{_html.escape(hess_path.name)}</span>'
            else:
                hess_html = f'<span style="color:#d32f2f;">{_html.escape(hess_path.name)} (missing)</span>'
            calc_print_mode_file_info.value = (
                f'<span style="color:#555;">Job:</span> <code>{_html.escape(out_path.name)}</code>'
                f' &nbsp;→&nbsp; '
                f'<span style="color:#555;">Hessian:</span> <code>{hess_html}</code>'
            )
        else:
            calc_print_mode_file_info.value = '<span style="color:#d32f2f;">Select a .out file.</span>'

    def _calc_show_print_mode_panel(show):
        if show:
            _calc_refresh_print_mode_info()
            calc_print_mode_panel.layout.display = 'block'
        else:
            calc_print_mode_panel.layout.display = 'none'
        _run_js(
            f"setTimeout(function(){{"
            f" if(window['{calc_resize_mol_fn}']) window['{calc_resize_mol_fn}']();"
            f"}}, 80);"
        )

    def _calc_mo_plot_parse_indices(raw):
        tokens = [tok.strip() for tok in re.split(r'[\s,]+', str(raw or '').strip()) if tok.strip()]
        if not tokens:
            return None, 'Please enter at least one MO index (e.g. 45 46 47).'
        indices = []
        for token in tokens:
            try:
                mo_idx = int(token)
            except ValueError:
                return None, f'Invalid MO index: {token}'
            if mo_idx < 0:
                return None, 'MO indices must be non-negative integers.'
            indices.append(mo_idx)
        return indices, ''

    def _calc_mo_project_dir(selected_path):
        if selected_path is None:
            return None
        try:
            calc_root = _calc_dir().resolve()
        except Exception:
            calc_root = _calc_dir()
        for parent in selected_path.parents:
            if (parent / 'CONTROL.txt').exists() or (parent / 'DELFIN_Data.json').exists():
                return parent
            try:
                if parent.resolve() == calc_root:
                    break
            except Exception:
                if parent == calc_root:
                    break
        return selected_path.parent

    def _calc_mo_find_gbw(selected_path):
        if selected_path is None or not selected_path.exists():
            return None

        base_dir = selected_path.parent
        stem = selected_path.stem
        stem_candidates = [stem]
        out_match = re.match(r'(?i)^output(\d*)$', stem)
        if out_match:
            mapped = f"input{out_match.group(1)}"
            if mapped not in stem_candidates:
                stem_candidates.append(mapped)
        inp_match = re.match(r'(?i)^input(\d*)$', stem)
        if inp_match:
            mapped = f"output{inp_match.group(1)}"
            if mapped not in stem_candidates:
                stem_candidates.append(mapped)

        candidates = []

        def _push(path):
            if path not in candidates:
                candidates.append(path)

        for st in stem_candidates:
            _push(base_dir / f'{st}.gbw')
        _push(base_dir / 'initial.gbw')
        _push(base_dir / 'ESD' / 'S0.gbw')
        try:
            for gbw in sorted(base_dir.glob('*.gbw'), key=lambda p: p.name.lower()):
                _push(gbw)
        except Exception:
            pass

        parent_dir = base_dir.parent
        if parent_dir != base_dir:
            for st in stem_candidates:
                _push(parent_dir / f'{st}.gbw')
            _push(parent_dir / 'initial.gbw')
            _push(parent_dir / 'ESD' / 'S0.gbw')
            try:
                for gbw in sorted(parent_dir.glob('*.gbw'), key=lambda p: p.name.lower()):
                    _push(gbw)
            except Exception:
                pass

        for path in candidates:
            if path.exists() and path.is_file():
                return path
        return None

    def _calc_mo_detect_spin(selected_path):
        scheme = None
        homo_index = None
        if selected_path is None or not selected_path.exists():
            return 'RKS', homo_index

        # Primary detection from selected .out (requested): HFTyp/Hartree-Fock type
        try:
            with selected_path.open('rb') as handle:
                head = handle.read(1200000).decode('utf-8', errors='ignore')
            hftyp_token = None
            m_hftyp = re.search(
                r'(?im)^\s*Hartree-?Fock\s+type.*?\b(UHF|RHF|ROHF|UKS|RKS)\b',
                head,
            )
            if m_hftyp:
                hftyp_token = str(m_hftyp.group(1)).upper().strip()
            if not hftyp_token:
                m_short = re.search(r'(?im)^\s*HFTyp\s*\.{2,}\s*([A-Za-z0-9_-]+)\b', head)
                if m_short:
                    hftyp_token = str(m_short.group(1)).upper().strip()
            if hftyp_token:
                scheme = 'UKS' if hftyp_token.startswith('U') else 'RKS'
        except Exception:
            head = ''

        project_dir = _calc_mo_project_dir(selected_path)
        json_candidates = []
        if project_dir is not None:
            json_candidates.append(project_dir / 'DELFIN_Data.json')
        for parent in selected_path.parents:
            candidate = parent / 'DELFIN_Data.json'
            if candidate not in json_candidates:
                json_candidates.append(candidate)
            if project_dir is not None and parent == project_dir:
                break

        for json_path in json_candidates:
            if not json_path.exists():
                continue
            try:
                data = json.loads(json_path.read_text(errors='ignore'))
            except Exception:
                continue
            orbitals = ((data.get('ground_state_S0') or {}).get('orbitals') or {})
            raw_scheme = str(orbitals.get('occupation_scheme') or '').upper().strip()
            if scheme is None:
                if raw_scheme in ('UKS', 'UHF'):
                    scheme = 'UKS'
                elif raw_scheme:
                    scheme = 'RKS'
            try:
                raw_homo = orbitals.get('homo_index')
                if raw_homo is not None and homo_index is None:
                    homo_index = int(raw_homo)
            except Exception:
                pass
            if scheme is not None and homo_index is not None:
                break

        # Fallback: parse orbitals directly from selected ORCA output.
        # This is robust for OCCUPIER directories without DELFIN_Data.json.
        try:
            from delfin.reporting.delfin_collector import parse_orbitals as _parse_orbitals

            parsed = _parse_orbitals(selected_path)
            if parsed:
                parsed_scheme = str(parsed.get('occupation_scheme') or '').upper().strip()
                if scheme is None:
                    if parsed_scheme in ('UKS', 'UHF'):
                        scheme = 'UKS'
                    elif parsed_scheme:
                        scheme = 'RKS'
                raw_homo = parsed.get('homo_index')
                if raw_homo is not None and homo_index is None:
                    homo_index = int(raw_homo)
        except Exception:
            pass

        try:
            tail = ''
            with selected_path.open('rb') as handle:
                handle.seek(0, 2)
                size = handle.tell()
                handle.seek(max(0, size - 300000))
                tail = handle.read().decode('utf-8', errors='ignore').upper()
            if scheme is None and (re.search(r'\b(UKS|UHF)\b', head.upper()) or re.search(r'\b(UKS|UHF)\b', tail)):
                scheme = 'UKS'
        except Exception:
            pass

        return (scheme or 'RKS'), homo_index

    def _calc_refresh_mo_plot_info():
        selected_path = _calc_get_selected_path()
        if selected_path is None or selected_path.suffix.lower() != '.out':
            calc_mo_plot_file_info.value = '<span style="color:#d32f2f;">Select a .out file.</span>'
            calc_mo_plot_spin_box.layout.display = 'none'
            return

        gbw_path = _calc_mo_find_gbw(selected_path)
        scheme, homo_index = _calc_mo_detect_spin(selected_path)
        scheme = 'UKS' if scheme == 'UKS' else 'RKS'
        calc_mo_plot_spin_box.layout.display = 'flex' if scheme == 'UKS' else 'none'

        if gbw_path is not None and gbw_path.exists():
            gbw_html = f'<span style="color:#2e7d32;">{_html.escape(str(gbw_path))}</span>'
        else:
            gbw_html = '<span style="color:#d32f2f;">not found</span>'
        homo_text = str(homo_index) if homo_index is not None else 'n/a'
        calc_mo_plot_file_info.value = (
            f'<span style="color:#555;">GBW:</span> {gbw_html}'
            f' &nbsp;|&nbsp; '
            f'<span style="color:#555;">Type:</span> <code>{_html.escape(scheme)}</code>'
            f' &nbsp;|&nbsp; '
            f'<span style="color:#555;">HOMO:</span> <code>{_html.escape(homo_text)}</code>'
        )

    def _calc_update_mo_spin_button_styles(_change=None):
        calc_mo_plot_alpha_cb.button_style = 'info' if calc_mo_plot_alpha_cb.value else ''
        calc_mo_plot_beta_cb.button_style = 'info' if calc_mo_plot_beta_cb.value else ''

    def _calc_show_mo_plot_panel(show):
        if show:
            _calc_refresh_mo_plot_info()
            _calc_update_mo_spin_button_styles()
            calc_mo_plot_panel.layout.display = 'block'
        else:
            calc_mo_plot_panel.layout.display = 'none'
        _run_js(
            f"setTimeout(function(){{"
            f" if(window['{calc_resize_mol_fn}']) window['{calc_resize_mol_fn}']();"
            f"}}, 80);"
        )

    def _calc_show_nmr_panel(show):
        if show:
            calc_nmr_status.value = ''
            calc_nmr_panel.layout.display = 'block'
        else:
            calc_nmr_panel.layout.display = 'none'

    def _calc_show_censo_nmr_panel(show):
        if show:
            calc_censo_nmr_status.value = ''
            calc_censo_nmr_panel.layout.display = 'block'
        else:
            calc_censo_nmr_panel.layout.display = 'none'

    def _calc_run_print_nmr():
        """Generate 1H NMR spectrum PNG from selected ORCA .out file."""
        selected_path = _calc_get_selected_path()
        if selected_path is None or not selected_path.exists():
            calc_content_area.value = '<span style="color:#d32f2f;">Select a .out file first.</span>'
            return
        if selected_path.suffix.lower() != '.out':
            calc_content_area.value = '<span style="color:#d32f2f;">Selected file is not a .out file.</span>'
            return

        calc_content_area.value = '<span style="color:#1565c0;">Generating NMR spectrum...</span>'

        def _run():
            try:
                from delfin.nmr_spectrum import parse_nmr_orca
                from delfin.reporting.nmr_report import (
                    create_nmr_report,
                    print_assignment_table,
                )

                result = parse_nmr_orca(selected_path)
                if not result.h_shieldings:
                    calc_content_area.value = (
                        '<span style="color:#d32f2f;">No hydrogen shieldings found '
                        'in this output file.</span>'
                    )
                    return

                out_png = selected_path.with_name(f'NMR_{selected_path.stem}.png')
                title = selected_path.stem.replace('_', ' ')
                create_nmr_report(
                    result,
                    output_png=out_png,
                    title=f'Calculated 1H NMR: {title}',
                )

                table_text = print_assignment_table(result)
                table_html = _html.escape(table_text).replace('\n', '<br>')

                calc_content_area.value = (
                    f'<span style="color:#2e7d32;">NMR spectrum saved: '
                    f'<code>{out_png.name}</code></span><br><br>'
                    f'<img src="files/{out_png.relative_to(selected_path.parent)}" '
                    f'style="max-width:100%;" /><br><br>'
                    f'<pre style="font-size:11px;">{table_html}</pre>'
                )
            except Exception as exc:
                import traceback as _tb
                calc_content_area.value = (
                    f'<span style="color:#d32f2f;">NMR generation failed: '
                    f'{_html.escape(str(exc))}</span>'
                    f'<pre style="font-size:10px;">{_html.escape(_tb.format_exc())}</pre>'
                )

        import threading
        threading.Thread(target=_run, daemon=True).start()

    def _calc_on_nmr_submit(_button):
        """Create ORCA NMR input from selected .xyz and submit SLURM job."""
        selected_path = _calc_selected_item_path()
        if selected_path is None or not selected_path.exists():
            calc_nmr_status.value = '<span style="color:#d32f2f;">Select a .xyz file first.</span>'
            return
        if selected_path.suffix.lower() != '.xyz':
            calc_nmr_status.value = '<span style="color:#d32f2f;">Selected file is not a .xyz file.</span>'
            return

        # Read xyz coordinates (skip first 2 lines: atom count + comment)
        try:
            xyz_lines = selected_path.read_text().splitlines()
            coord_lines = xyz_lines[2:]  # skip atom count + comment
            # Filter empty lines at the end
            while coord_lines and not coord_lines[-1].strip():
                coord_lines.pop()
            if not coord_lines:
                calc_nmr_status.value = '<span style="color:#d32f2f;">XYZ file is empty.</span>'
                return
        except Exception as exc:
            calc_nmr_status.value = (
                f'<span style="color:#d32f2f;">Failed to read XYZ: {_html.escape(str(exc))}</span>'
            )
            return

        pal = max(1, int(calc_nmr_pal.value or 12))
        maxcore = max(100, int(calc_nmr_maxcore.value or 3000))
        time_limit = calc_nmr_time.value.strip() or '03:00:00'
        solvent = str(calc_nmr_solvent.value or 'chloroform').strip() or 'chloroform'

        # Create NMR subfolder
        nmr_dir = _calc_next_available_dir(selected_path.parent, f'{selected_path.stem}_NMR')
        job_name = nmr_dir.name
        try:
            nmr_dir.mkdir(parents=True, exist_ok=False)
        except Exception as exc:
            calc_nmr_status.value = (
                f'<span style="color:#d32f2f;">Failed to create folder: {_html.escape(str(exc))}</span>'
            )
            return

        # Generate .inp content
        inp_content = build_calc_nmr_input(
            coord_lines,
            pal=pal,
            maxcore=maxcore,
            solvent=solvent,
        )

        inp_file = f'{job_name}.inp'
        try:
            (nmr_dir / inp_file).write_text(inp_content)
        except Exception as exc:
            calc_nmr_status.value = (
                f'<span style="color:#d32f2f;">Failed to write .inp: {_html.escape(str(exc))}</span>'
            )
            return

        # Submit
        try:
            result = ctx.backend.submit_orca(
                job_dir=nmr_dir,
                job_name=job_name,
                inp_file=inp_file,
                time_limit=time_limit,
                pal=pal,
                maxcore=maxcore,
            )
        except Exception as exc:
            calc_nmr_status.value = (
                f'<span style="color:#d32f2f;">Submit error: {_html.escape(str(exc))}</span>'
            )
            return

        if result.returncode == 0:
            job_id = result.stdout.strip().split()[-1] if result.stdout.strip() else '(unknown)'
            calc_nmr_status.value = (
                f'<span style="color:#2e7d32;">Submitted <code>{_html.escape(job_name)}</code>'
                f' (ID: {_html.escape(job_id)}) with <code>CPCM({_html.escape(solvent)})</code>'
                f' &mdash; .inp in {_html.escape(nmr_dir.name)}/</span>'
            )
        else:
            msg = result.stderr or result.stdout or 'Unknown error'
            calc_nmr_status.value = f'<span style="color:#d32f2f;">{_html.escape(msg)}</span>'

        _calc_show_nmr_panel(False)

    def _calc_on_censo_nmr_submit(_button):
        """Submit a CREST + CENSO + ANMR ensemble NMR workflow for the selected XYZ."""
        selected_path = _calc_selected_item_path()
        if selected_path is None or not selected_path.exists():
            calc_censo_nmr_status.value = '<span style="color:#d32f2f;">Select a .xyz file first.</span>'
            return
        if selected_path.suffix.lower() != '.xyz':
            calc_censo_nmr_status.value = '<span style="color:#d32f2f;">Selected file is not a .xyz file.</span>'
            return

        pal = max(1, int(calc_censo_nmr_pal.value or 12))
        maxcore = max(100, int(calc_censo_nmr_maxcore.value or 3000))
        time_limit = calc_censo_nmr_time.value.strip() or '12:00:00'
        solvent = str(calc_censo_nmr_solvent.value or 'chcl3').strip().lower() or 'chcl3'
        charge = int(calc_censo_nmr_charge.value or 0)
        multiplicity = max(1, int(calc_censo_nmr_multiplicity.value or 1))
        mhz = max(50, int(calc_censo_nmr_mhz.value or 400))

        job_dir = _calc_next_available_dir(selected_path.parent, f'{selected_path.stem}_CENSO_ANMR')
        job_name = job_dir.name
        try:
            job_dir.mkdir(parents=True, exist_ok=False)
            source_copy = job_dir / selected_path.name
            shutil.copy2(selected_path, source_copy)
        except Exception as exc:
            calc_censo_nmr_status.value = (
                f'<span style="color:#d32f2f;">Failed to prepare workflow folder: {_html.escape(str(exc))}</span>'
            )
            return

        try:
            result = ctx.backend.submit_delfin(
                job_dir=job_dir,
                job_name=job_name,
                mode='censo_anmr',
                time_limit=time_limit,
                pal=pal,
                maxcore=maxcore,
                extra_env={
                    'DELFIN_XYZ_FILE': str(source_copy),
                    'DELFIN_WORKFLOW_LABEL': job_name,
                    'DELFIN_CENSO_NMR_SOLVENT': solvent,
                    'DELFIN_CENSO_NMR_CHARGE': str(charge),
                    'DELFIN_CENSO_NMR_MULTIPLICITY': str(multiplicity),
                    'DELFIN_CENSO_NMR_MHZ': str(mhz),
                },
            )
        except Exception as exc:
            calc_censo_nmr_status.value = (
                f'<span style="color:#d32f2f;">Submit error: {_html.escape(str(exc))}</span>'
            )
            return

        if result.returncode == 0:
            job_id = result.stdout.strip().split()[-1] if result.stdout.strip() else '(unknown)'
            calc_censo_nmr_status.value = (
                f'<span style="color:#2e7d32;">Submitted <code>{_html.escape(job_name)}</code>'
                f' (ID: {_html.escape(job_id)}) for <code>{_html.escape(solvent)}</code>'
                f' ensemble 1H NMR in <code>{_html.escape(job_dir.name)}</code>.</span>'
            )
            _calc_show_censo_nmr_panel(False)
            return

        msg = result.stderr or result.stdout or 'Unknown error'
        calc_censo_nmr_status.value = f'<span style="color:#d32f2f;">{_html.escape(msg)}</span>'

    def _calc_render_traj_plot():
        selected_path = _calc_get_selected_path()
        if selected_path is None or not selected_path.exists():
            calc_download_status.value = (
                '<span style="color:#d32f2f;">Select a .final.interp file first.</span>'
            )
            return
        try:
            raw_text = state.get('file_content') or selected_path.read_text(errors='ignore')
            output_path_hartree = selected_path.with_name(f'{selected_path.stem}_hartree.png')
            output_path_kcal = selected_path.with_name(f'{selected_path.stem}_kcal_mol.png')
            csv_output_path_point = selected_path.with_name(f'{selected_path.stem}_point.csv')
            csv_output_path_comma = selected_path.with_name(f'{selected_path.stem}_comma.csv')
            save_neb_trajectory_plot_png(
                raw_text,
                output_path_hartree,
                title=selected_path.name,
                energy_unit='hartree',
            )
            save_neb_trajectory_plot_png(
                raw_text,
                output_path_kcal,
                title=selected_path.name,
                energy_unit='kcal/mol',
            )
            save_neb_trajectory_csv(raw_text, csv_output_path_point, decimal='.')
            save_neb_trajectory_csv(raw_text, csv_output_path_comma, decimal=',')
            calc_download_status.value = (
                '<span style="color:#2e7d32;">Saved trajectory files:</span> '
                f'<code>{_html.escape(output_path_hartree.name)}</code>'
                ' &nbsp;and&nbsp; '
                f'<code>{_html.escape(output_path_kcal.name)}</code>'
                ' &nbsp;and&nbsp; '
                f'<code>{_html.escape(csv_output_path_point.name)}</code>'
                ' &nbsp;and&nbsp; '
                f'<code>{_html.escape(csv_output_path_comma.name)}</code>'
                ' <span style="color:#555;">(CSV uses ; as separator; both contain Eh and kcal/mol)</span>'
            )
        except Exception as exc:
            calc_download_status.value = f'<span style="color:#d32f2f;">{_html.escape(str(exc))}</span>'

    def calc_on_print_mode_plot(_button):
        modes, err = _calc_print_mode_parse_modes(calc_print_mode_input.value)
        if err:
            calc_print_mode_status.value = f'<span style="color:#d32f2f;">{_html.escape(err)}</span>'
            return

        out_path, hess_path, path_err = _calc_print_mode_target_paths()
        _calc_refresh_print_mode_info()
        if path_err:
            calc_print_mode_status.value = (
                f'<span style="color:#d32f2f;">{_html.escape(path_err)}</span>'
            )
            return

        if shutil.which('orca_pltvib') is None:
            calc_print_mode_status.value = (
                '<span style="color:#d32f2f;">orca_pltvib not found in PATH.</span>'
            )
            return

        cmd = ['orca_pltvib', hess_path.name, *modes]
        cmd_text = _html.escape(' '.join(cmd))
        workdir_text = _html.escape(str(hess_path.parent))
        calc_print_mode_status.value = (
            f'<span style="color:#1976d2;">Starting: <code>{cmd_text}</code> '
            f'in <code>{workdir_text}</code> ...</span>'
        )

        def _run_print_mode():
            try:
                result = subprocess.run(
                    cmd,
                    cwd=str(hess_path.parent),
                    capture_output=True,
                    text=True,
                    check=False,
                )
                if result.returncode == 0:
                    calc_print_mode_status.value = (
                        f'<span style="color:#2e7d32;">Finished: <code>{cmd_text}</code> '
                        f'for <code>{_html.escape(out_path.name)}</code>.</span>'
                    )
                else:
                    err = (result.stderr or result.stdout or '').strip()
                    if len(err) > 240:
                        err = err[:240] + '...'
                    calc_print_mode_status.value = (
                        '<span style="color:#d32f2f;">'
                        f'orca_pltvib failed (rc={result.returncode}): '
                        f'{_html.escape(err or "no error output")}'
                        '</span>'
                    )
            except Exception as exc:
                calc_print_mode_status.value = (
                    f'<span style="color:#d32f2f;">Failed: {_html.escape(str(exc))}</span>'
                )

        threading.Thread(target=_run_print_mode, daemon=True).start()
        calc_options_dropdown.value = '(Options)'

    def calc_on_mo_plot(_button):
        mo_indices, err = _calc_mo_plot_parse_indices(calc_mo_plot_input.value)
        if err:
            calc_mo_plot_status.value = f'<span style="color:#d32f2f;">{_html.escape(err)}</span>'
            return

        selected_path = _calc_get_selected_path()
        _calc_refresh_mo_plot_info()
        if selected_path is None or selected_path.suffix.lower() != '.out':
            calc_mo_plot_status.value = '<span style="color:#d32f2f;">Select a .out file first.</span>'
            return

        gbw_path = _calc_mo_find_gbw(selected_path)
        if gbw_path is None:
            calc_mo_plot_status.value = (
                '<span style="color:#d32f2f;">No suitable .gbw file found for selected job.</span>'
            )
            return

        orca_plot_exe = shutil.which('orca_plot') or _orca_plot_binary()
        if not orca_plot_exe:
            calc_mo_plot_status.value = (
                '<span style="color:#d32f2f;">orca_plot not found in PATH or ORCA env variables.</span>'
            )
            return

        scheme, _homo_index = _calc_mo_detect_spin(selected_path)
        scheme = 'UKS' if scheme == 'UKS' else 'RKS'
        if scheme == 'UKS':
            spin_choices = []
            if calc_mo_plot_alpha_cb.value:
                spin_choices.append(('alpha', 0))
            if calc_mo_plot_beta_cb.value:
                spin_choices.append(('beta', 1))
            if not spin_choices:
                calc_mo_plot_status.value = (
                    '<span style="color:#d32f2f;">Select Alpha and/or Beta for UKS jobs.</span>'
                )
                return
        else:
            spin_choices = [('alpha', 0)]

        jobs = [(mo_idx, spin_name, spin_op) for mo_idx in mo_indices for spin_name, spin_op in spin_choices]
        gbw_dir = gbw_path.parent
        project_dir = _calc_mo_project_dir(selected_path) or gbw_dir
        calc_mo_plot_btn.disabled = True
        calc_mo_plot_status.value = (
            f'<span style="color:#1976d2;">Generating {len(jobs)} cube job(s) with '
            f'<code>{_html.escape(str(orca_plot_exe))}</code> ...</span>'
        )

        def _run_mo_plot():
            generated = []
            failures = []
            created_temp_files = []
            env = os.environ.copy()
            env['ORCA_SCRDIR'] = str(gbw_dir)
            env['ORCA_TMPDIR'] = str(gbw_dir)
            env.setdefault('TMPDIR', str(gbw_dir))

            def _collect_cube_files():
                found = []
                try:
                    found.extend(gbw_dir.glob('*.cube'))
                except Exception:
                    pass
                try:
                    found.extend(gbw_dir.glob('*.cub'))
                except Exception:
                    pass
                uniq = []
                seen = set()
                for p in found:
                    key = str(p.resolve())
                    if key in seen:
                        continue
                    seen.add(key)
                    uniq.append(p)
                return sorted(uniq, key=lambda p: p.stat().st_mtime_ns)

            try:
                dens_candidates = [
                    project_dir / f'{gbw_path.stem}.densitiesinfo',
                    gbw_dir / f'{gbw_path.stem}.densitiesinfo',
                    project_dir / 'initial.densitiesinfo',
                    gbw_dir / 'initial.densitiesinfo',
                ]
                dens_src = next((p for p in dens_candidates if p.exists()), None)
                if dens_src is not None:
                    for dens_dest in {gbw_dir / dens_src.name, gbw_dir / f'{gbw_path.stem}.densitiesinfo'}:
                        if dens_dest.exists():
                            continue
                        try:
                            shutil.copy2(dens_src, dens_dest)
                            created_temp_files.append(dens_dest)
                        except Exception:
                            pass

                total = len(jobs)
                for idx, (mo_idx, spin_name, spin_op) in enumerate(jobs, start=1):
                    display_spin = spin_name if scheme == 'UKS' else 'RKS'
                    calc_mo_plot_status.value = (
                        '<span style="color:#1976d2;">'
                        f'Generating MO {mo_idx} ({_html.escape(display_spin)}) [{idx}/{total}] ...'
                        '</span>'
                    )

                    orca_input = f"1\n1\n2\n{mo_idx}\n3\n{spin_op}\n4\n70\n11\n12\n"
                    cubes_before = _collect_cube_files()
                    before_state = {}
                    for prev_cube in cubes_before:
                        try:
                            st_prev = prev_cube.stat()
                        except Exception:
                            continue
                        before_state[str(prev_cube.resolve())] = (
                            st_prev.st_mtime_ns,
                            st_prev.st_ctime_ns,
                            st_prev.st_size,
                        )
                    try:
                        result = subprocess.run(
                            [orca_plot_exe, gbw_path.name, '-i'],
                            input=orca_input.encode(),
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE,
                            cwd=str(gbw_dir),
                            env=env,
                            timeout=120,
                            check=False,
                        )
                    except subprocess.TimeoutExpired:
                        failures.append(f'MO {mo_idx} ({display_spin}): timeout after 120s')
                        continue
                    except Exception as exc:
                        failures.append(f'MO {mo_idx} ({display_spin}): {exc}')
                        continue

                    if result.returncode != 0:
                        stderr_raw = (result.stderr or result.stdout or b'').decode(errors='ignore').strip()
                        retry_path = None
                        match = re.search(r"Filename:\s*(\S+\.densitiesinfo)", stderr_raw)
                        if match and dens_src is not None:
                            retry_path = Path(match.group(1)).expanduser()
                            if not retry_path.is_absolute():
                                retry_path = gbw_dir / retry_path
                        if retry_path and dens_src is not None:
                            try:
                                retry_path.parent.mkdir(parents=True, exist_ok=True)
                                if not retry_path.exists():
                                    shutil.copy2(dens_src, retry_path)
                                    created_temp_files.append(retry_path)
                                result = subprocess.run(
                                    [orca_plot_exe, gbw_path.name, '-i'],
                                    input=orca_input.encode(),
                                    stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE,
                                    cwd=str(gbw_dir),
                                    env=env,
                                    timeout=120,
                                    check=False,
                                )
                                stderr_raw = (result.stderr or result.stdout or b'').decode(errors='ignore').strip()
                            except Exception:
                                pass
                    if result.returncode != 0:
                        stderr_text = (result.stderr or result.stdout or b'').decode(errors='ignore').strip()
                        stderr_text = stderr_text.replace('\n', ' ')
                        if len(stderr_text) > 220:
                            stderr_text = stderr_text[:220] + '...'
                        failures.append(
                            f'MO {mo_idx} ({display_spin}): rc={result.returncode} '
                            f'{stderr_text or "no error output"}'
                        )
                        continue

                    changed_nonempty = []
                    changed_empty = []
                    for _poll in range(12):
                        changed_nonempty = []
                        changed_empty = []
                        cube_files = _collect_cube_files()
                        for cube_file in cube_files:
                            try:
                                st_now = cube_file.stat()
                            except Exception:
                                continue
                            key = str(cube_file.resolve())
                            prev = before_state.get(key)
                            changed = (
                                prev is None
                                or st_now.st_mtime_ns > prev[0]
                                or st_now.st_ctime_ns > prev[1]
                                or st_now.st_size != prev[2]
                            )
                            if not changed:
                                continue
                            entry = (st_now.st_mtime_ns, cube_file, st_now.st_size)
                            if st_now.st_size > 0:
                                changed_nonempty.append(entry)
                            else:
                                changed_empty.append(entry)
                        if changed_nonempty:
                            break
                        if _poll < 11:
                            time.sleep(0.2)

                    if not changed_nonempty and not changed_empty:
                        failures.append(f'MO {mo_idx} ({display_spin}): no cube file generated')
                        continue
                    if not changed_nonempty:
                        empty_names = ', '.join(item[1].name for item in sorted(changed_empty)[-2:])
                        failures.append(
                            f'MO {mo_idx} ({display_spin}): only empty cube file(s) generated ({empty_names})'
                        )
                        continue

                    changed_nonempty.sort(key=lambda item: item[0])
                    cube_source = changed_nonempty[-1][1]
                    generated.append(cube_source)

                if generated:
                    shown = ', '.join(path.name for path in generated[:6])
                    if len(generated) > 6:
                        shown += f', ... (+{len(generated) - 6})'
                    if failures:
                        failure_text = '; '.join(failures[:2])
                        if len(failures) > 2:
                            failure_text += f' ... (+{len(failures) - 2} more)'
                        calc_mo_plot_status.value = (
                            '<span style="color:#ef6c00;">'
                            f'Generated {len(generated)} cube(s): {_html.escape(shown)}'
                            f' &nbsp;|&nbsp; Issues: {_html.escape(failure_text)}'
                            '</span>'
                        )
                    else:
                        calc_mo_plot_status.value = (
                            '<span style="color:#2e7d32;">'
                            f'Generated {len(generated)} cube(s): {_html.escape(shown)}'
                            '</span>'
                        )
                else:
                    failure_text = '; '.join(failures[:2]) if failures else 'No cube files were generated.'
                    if failures and len(failures) > 2:
                        failure_text += f' ... (+{len(failures) - 2} more)'
                    calc_mo_plot_status.value = (
                        f'<span style="color:#d32f2f;">{_html.escape(failure_text)}</span>'
                    )
            finally:
                for temp_path in created_temp_files:
                    try:
                        if temp_path.exists():
                            temp_path.unlink()
                    except Exception:
                        pass
                calc_mo_plot_btn.disabled = False

        threading.Thread(target=_run_mo_plot, daemon=True).start()
        calc_options_dropdown.value = '(Options)'

    def calc_on_xyz_batch_refresh(_button=None):
        _calc_refresh_xyz_batch_selector()
        calc_xyz_batch_status.value = '<span style="color:#555;">Selection refreshed.</span>'

    def _calc_toggle_xyz_batch_mark_values(values):
        selected_values = []
        seen = set()
        for value in values or ():
            item = str(value).strip()
            if not item or item in seen:
                continue
            seen.add(item)
            selected_values.append(item)
        if not selected_values:
            return 0, ''
        marked = set(state.get('xyz_batch_marked_paths') or [])
        if all(value in marked for value in selected_values):
            for value in selected_values:
                marked.discard(value)
            action = 'Removed'
        else:
            for value in selected_values:
                marked.add(value)
            action = 'Marked'
        state['xyz_batch_marked_paths'] = sorted(marked)
        _calc_refresh_xyz_batch_selector()
        return len(selected_values), action

    def calc_on_xyz_batch_mark(_button=None):
        selected_values = [str(value).strip() for value in (calc_xyz_batch_select.value or ()) if str(value).strip()]
        if not selected_values:
            calc_xyz_batch_status.value = (
                '<span style="color:#d32f2f;">Select .xyz files or folders first, then click Select.</span>'
            )
            return
        count, action = _calc_toggle_xyz_batch_mark_values(selected_values)
        calc_xyz_batch_status.value = (
            f'<span style="color:#2e7d32;">{_html.escape(action)} {count} item(s).</span>'
        )

    def calc_on_xyz_batch_toggle(change):
        raw = (change['new'] or '').strip()
        calc_xyz_batch_toggle_input.value = ''
        if not raw:
            return
        toggle_value = raw
        if raw.startswith('idx:'):
            try:
                idx = int(raw.split(':', 1)[1])
            except Exception:
                idx = -1
            candidates = _calc_collect_xyz_batch_labels(_calc_xyz_batch_dir())
            if 0 <= idx < len(candidates):
                toggle_value = str(candidates[idx][1]).strip()
            else:
                toggle_value = ''
        count, action = _calc_toggle_xyz_batch_mark_values([toggle_value] if toggle_value else [])
        if count <= 0:
            return
        calc_xyz_batch_status.value = (
            f'<span style="color:#2e7d32;">{_html.escape(action)} 1 item.</span>'
        )

    def calc_on_xyz_batch_select_all(_button=None):
        candidates = _calc_collect_xyz_batch_labels(_calc_xyz_batch_dir())
        visible_values = [value for _label, value in candidates]
        if not visible_values:
            calc_xyz_batch_status.value = (
                '<span style="color:#d32f2f;">No .xyz files or folders in this directory.</span>'
            )
            return
        marked = set(state.get('xyz_batch_marked_paths') or [])
        marked.update(visible_values)
        state['xyz_batch_marked_paths'] = sorted(marked)
        calc_xyz_batch_select.value = tuple(visible_values)
        _calc_refresh_xyz_batch_selector()
        calc_xyz_batch_status.value = (
            f'<span style="color:#2e7d32;">Marked all {len(visible_values)} visible item(s) in this directory.</span>'
        )

    def calc_on_xyz_batch_up(_button=None):
        current_dir = _calc_xyz_batch_dir()
        root_dir = _calc_dir()
        if current_dir == root_dir:
            calc_xyz_batch_status.value = '<span style="color:#555;">Already at calc root.</span>'
            return
        if not _calc_set_xyz_batch_dir(current_dir.parent):
            calc_xyz_batch_status.value = '<span style="color:#d32f2f;">Could not open parent folder.</span>'
            return
        _calc_refresh_xyz_batch_selector()
        calc_xyz_batch_status.value = '<span style="color:#555;">Moved to parent folder.</span>'

    def calc_on_xyz_batch_root(_button=None):
        if not _calc_set_xyz_batch_dir(_calc_dir()):
            calc_xyz_batch_status.value = '<span style="color:#d32f2f;">Could not open calc root.</span>'
            return
        _calc_refresh_xyz_batch_selector()
        calc_xyz_batch_status.value = '<span style="color:#555;">Moved to calc root.</span>'

    def calc_on_xyz_batch_dblclick(change):
        raw = (change['new'] or '').strip()
        calc_xyz_batch_dblclick_input.value = ''
        if not raw or raw.startswith('('):
            return
        item_name = _calc_label_to_name(raw)
        if not item_name:
            return
        item_path = _calc_xyz_batch_dir() / item_name
        if not item_path.exists():
            calc_xyz_batch_status.value = (
                f'<span style="color:#d32f2f;">{_html.escape(item_name)} not found.</span>'
            )
            return
        if not item_path.is_dir():
            return
        if not _calc_set_xyz_batch_dir(item_path):
            calc_xyz_batch_status.value = (
                f'<span style="color:#d32f2f;">Could not open folder '
                f'<code>{_html.escape(item_name)}</code>.</span>'
            )
            return
        _calc_refresh_xyz_batch_selector()
        calc_xyz_batch_status.value = (
            f'<span style="color:#555;">Opened folder <code>{_html.escape(item_name)}</code>.</span>'
        )

    def _calc_xyz_batch_skipped_html(skipped):
        if not skipped:
            return ''
        details = '; '.join(_html.escape(s) for s in skipped[:3])
        more = '' if len(skipped) <= 3 else f' (+{len(skipped) - 3} more)'
        return f' <span style="color:#ef6c00;">Skipped {len(skipped)}: {details}{more}</span>'

    def _calc_xyz_batch_prepare_export():
        marked_values = [str(value).strip() for value in (state.get('xyz_batch_marked_paths') or []) if str(value).strip()]
        selected_values = [str(value).strip() for value in (calc_xyz_batch_select.value or ()) if str(value).strip()]
        item_values = marked_values or selected_values
        if not item_values:
            return None, [], 'No files/jobs selected. Use Haken or Select All.'

        requested = calc_xyz_batch_filename.value.strip() or _calc_xyz_batch_default_filename()
        batch_stem = re.sub(r'[^A-Za-z0-9._-]+', '_', requested).strip('._')
        if not batch_stem:
            batch_stem = _calc_xyz_batch_default_filename()
        if batch_stem.lower().endswith('.txt'):
            batch_stem = batch_stem[:-4].rstrip('._') or 'batch'

        entries = []
        skipped = []
        home_dir = Path.home().resolve()
        for rel_value in item_values:
            try:
                item_path = _calc_resolve_within_root(_calc_dir() / rel_value)
            except Exception:
                skipped.append(f'{rel_value}: invalid item path')
                continue
            item_name = item_path.name
            if not item_path.exists():
                skipped.append(f'{item_name}: item not found')
                continue

            xyz_path, xyz_err = _calc_xyz_batch_resolve_source(item_path)
            if xyz_path is None:
                skipped.append(xyz_err or f'{item_name}: no xyz source')
                continue

            try:
                raw_xyz = xyz_path.read_text(errors='ignore')
            except Exception as exc:
                skipped.append(f'{item_name}: read failed ({exc})')
                continue

            xyz_content = _calc_xyz_batch_clean_xyz(raw_xyz)
            if not xyz_content:
                skipped.append(f'{item_name}: empty xyz coordinates')
                continue

            base_name = item_path.name if item_path.is_dir() else item_path.stem
            safe_name = ''.join(c for c in str(base_name) if c.isalnum() or c in ('_', '-'))
            if not safe_name:
                skipped.append(f'{item_name}: invalid batch name')
                continue
            entry_name = f'{safe_name}_{batch_stem}' if batch_stem else safe_name
            try:
                xyz_source_path = xyz_path.resolve()
            except Exception:
                xyz_source_path = xyz_path
            try:
                path_field = str(xyz_source_path.relative_to(_calc_dir()))
            except Exception:
                try:
                    path_field = str(xyz_source_path.relative_to(home_dir))
                except Exception:
                    path_field = str(xyz_source_path)
            entries.append((entry_name, path_field, xyz_content))

        if not entries:
            return None, skipped, 'No valid XYZ entries found.'

        lines = []
        for idx, (entry_name, rel_value, xyz_content) in enumerate(entries):
            lines.append(f'{entry_name};source={rel_value};')
            lines.extend(xyz_content.splitlines())
            lines.append('*')
            if idx < len(entries) - 1:
                lines.append('')
        payload_text = '\n'.join(lines).rstrip() + '\n'

        safe_file = re.sub(r'[^A-Za-z0-9._-]+', '_', requested).strip('._')
        if not safe_file:
            safe_file = _calc_xyz_batch_default_filename()
        if not safe_file.lower().endswith('.txt'):
            safe_file += '.txt'

        export_data = {
            'payload_text': payload_text,
            'safe_file': safe_file,
            'entry_count': len(entries),
        }
        return export_data, skipped, ''

    def calc_on_xyz_batch_build(_button):
        export_data, skipped, error = _calc_xyz_batch_prepare_export()
        if error:
            calc_xyz_batch_status.value = (
                f'<span style="color:#d32f2f;">{_html.escape(error)}</span>'
                + _calc_xyz_batch_skipped_html(skipped)
            )
            return

        payload = export_data['payload_text'].encode('utf-8')
        safe_file = export_data['safe_file']
        entry_count = int(export_data['entry_count'])
        _calc_trigger_download(safe_file, payload, mime='text/plain;charset=utf-8')

        status = (
            f'<span style="color:#2e7d32;">Download started: '
            f'{_html.escape(safe_file)} ({entry_count} entries)</span>'
        )
        status += _calc_xyz_batch_skipped_html(skipped)
        calc_xyz_batch_status.value = status

    def calc_on_xyz_batch_copy(_button):
        export_data, skipped, error = _calc_xyz_batch_prepare_export()
        if error:
            calc_xyz_batch_status.value = (
                f'<span style="color:#d32f2f;">{_html.escape(error)}</span>'
                + _calc_xyz_batch_skipped_html(skipped)
            )
            return

        payload_text = export_data['payload_text']
        safe_file = export_data['safe_file']
        entry_count = int(export_data['entry_count'])
        _calc_copy_to_clipboard(payload_text, label='batch xyz export')
        status = (
            f'<span style="color:#2e7d32;">Copy prepared (paste with Cmd+V/Ctrl+V): '
            f'{_html.escape(safe_file)} ({entry_count} entries)</span>'
        )
        status += _calc_xyz_batch_skipped_html(skipped)
        calc_xyz_batch_status.value = status

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

    def _calc_clear_plain_match_marker():
        _run_js("""
        setTimeout(function() {
            const old = document.getElementById('calc-current-match');
            if (!old) return;
            const parent = old.parentNode;
            if (!parent) return;
            parent.replaceChild(document.createTextNode(old.textContent || ''), old);
            parent.normalize();
        }, 0);
        """)

    def _calc_focus_plain_match(start_pos, end_pos):
        js = r"""
        setTimeout(function() {
            const box = document.getElementById('calc-content-box');
            const el = document.getElementById('calc-content-text');
            if (!box || !el) return;

            function unwrapCurrent() {
                const old = document.getElementById('calc-current-match');
                if (!old) return;
                const parent = old.parentNode;
                if (!parent) return;
                parent.replaceChild(document.createTextNode(old.textContent || ''), old);
                parent.normalize();
            }

            function scrollByRatio(startIndex, textLength) {
                if (textLength <= 0) return;
                const ratio = Math.max(0, Math.min(1, startIndex / Math.max(1, textLength - 1)));
                const maxScroll = Math.max(0, box.scrollHeight - box.clientHeight);
                box.scrollTop = Math.floor(ratio * maxScroll);
            }

            unwrapCurrent();
            const text = el.textContent || '';
            const s = __START__;
            const e = __END__;
            if (s < 0 || e <= s || s >= text.length) return;

            function locate(root, index) {
                let remaining = index;
                const walker = document.createTreeWalker(root, NodeFilter.SHOW_TEXT, null);
                let node = walker.nextNode();
                while (node) {
                    const len = (node.nodeValue || '').length;
                    if (remaining <= len) {
                        return {node: node, offset: remaining};
                    }
                    remaining -= len;
                    node = walker.nextNode();
                }
                return null;
            }

            const endBound = Math.min(e, text.length);
            const startLoc = locate(el, s);
            const endLoc = locate(el, endBound);
            if (!startLoc || !endLoc) {
                scrollByRatio(s, text.length);
                return;
            }

            const range = document.createRange();
            range.setStart(startLoc.node, startLoc.offset);
            range.setEnd(endLoc.node, endLoc.offset);

            const mark = document.createElement('mark');
            mark.className = 'calc-match current';
            mark.id = 'calc-current-match';
            try {
                range.surroundContents(mark);
            } catch (_err) {
                scrollByRatio(s, text.length);
                return;
            }

            mark.scrollIntoView({block: 'center'});
        }, 0);
        """
        js = js.replace('__START__', str(int(start_pos))).replace('__END__', str(int(end_pos)))
        _run_js(js)

    def calc_update_view():
        if _is_archive_tab and state.get('table_panel_active', False):
            calc_mol_container.layout.display = 'none'
            calc_content_area.layout.display = 'none'
            calc_edit_area.layout.display = 'none'
            calc_content_label.layout.display = 'none'
            calc_content_toolbar.layout.display = 'none'
            calc_recalc_toolbar.layout.display = 'none'
            calc_xyz_workflow_toolbar.layout.display = 'none'
            calc_nmr_panel.layout.display = 'none'
            calc_censo_nmr_panel.layout.display = 'none'
            return
        show_mol = calc_view_toggle.value
        if show_mol:
            calc_mol_container.layout.display = 'block'
            calc_content_area.layout.display = 'none'
            calc_edit_area.layout.display = 'none'
            calc_content_label.layout.display = 'none'
            calc_content_toolbar.layout.display = 'none'
            calc_recalc_toolbar.layout.display = 'none'
            calc_xyz_workflow_toolbar.layout.display = 'none'
            calc_nmr_panel.layout.display = 'none'
            calc_censo_nmr_panel.layout.display = 'none'
        else:
            _calc_stop_xyz_playback(update_button=True)
            calc_mol_container.layout.display = 'none'
            calc_content_label.layout.display = 'block'
            if state['recalc_active']:
                calc_content_toolbar.layout.display = 'none'
                calc_content_area.layout.display = 'none'
                calc_edit_area.layout.display = 'block'
                calc_recalc_toolbar.layout.display = 'flex'
                calc_xyz_workflow_toolbar.layout.display = 'none'
            elif state['xyz_workflow_active']:
                calc_content_toolbar.layout.display = 'none'
                calc_content_area.layout.display = 'block'
                calc_edit_area.layout.display = 'none'
                calc_recalc_toolbar.layout.display = 'none'
                calc_xyz_workflow_toolbar.layout.display = 'flex'
            else:
                calc_content_toolbar.layout.display = 'flex'
                calc_content_area.layout.display = 'block'
                calc_edit_area.layout.display = 'none'
                calc_recalc_toolbar.layout.display = 'none'
                calc_xyz_workflow_toolbar.layout.display = 'none'
            calc_nmr_panel.layout.display = 'none'
            calc_censo_nmr_panel.layout.display = 'none'

    def calc_set_message(message):
        calc_content_area.value = (
            "<div style='height:100%; overflow-y:auto;"
            " border:1px solid #ddd; padding:6px; font-family:monospace;"
            " white-space:pre; background:#fafafa;'>"
            f"{_html.escape(message)}"
            "</div>"
        )

    _calc_path_syncing = [False]

    def calc_update_path_display():
        display_path = '/' + state['current_path'] if state['current_path'] else '/'
        _calc_path_syncing[0] = True
        try:
            calc_path_input.value = display_path
        finally:
            _calc_path_syncing[0] = False
        calc_back_btn.disabled = (state['current_path'] == '')

    def calc_update_report_btn():
        calc_report_btn.disabled = not bool(calc_collect_report_targets(limit=1))

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
                    if (mark) { mark.scrollIntoView({block: 'center'}); }
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
            _calc_stop_xyz_playback(update_button=True)
            saved = {}
            widgets_to_hide = [
                ('calc_mol_header', calc_mol_header),
                ('calc_xyz_frame_label', calc_xyz_frame_label),
                ('calc_xyz_controls', calc_xyz_controls),
                ('calc_coord_controls', calc_coord_controls),
                ('calc_mol_view_row', calc_mol_view_row),
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
            ('calc_mol_header', calc_mol_header),
            ('calc_xyz_frame_label', calc_xyz_frame_label),
            ('calc_xyz_controls', calc_xyz_controls),
            ('calc_coord_controls', calc_coord_controls),
            ('calc_mol_view_row', calc_mol_view_row),
            ('calc_mol_viewer', calc_mol_viewer),
        ]
        for key, widget in widgets_to_restore:
            widget.layout.display = saved.get(key, '')
        state['rmsd_saved_display'] = {}
        state['rmsd_mode_active'] = False
        _calc_update_xyz_traj_control_state()

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
                var rightPanel = scopeRoot ? scopeRoot.querySelector('.calc-right') : null;
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
                                if (!style || style.display === 'none' || style.visibility === 'hidden') {{
                                    continue;
                                }}
                                var childRect = child.getBoundingClientRect();
                                if (childRect.height > 0) reservedBelow += childRect.height;
                            }}
                        }}
                        var availH = rightRect.bottom - mvRect.top - reservedBelow - 10;
                        var row = mv.closest('.calc-mol-view-row');
                        var rowRect = row ? row.getBoundingClientRect() : rightRect;
                        var tray = scopeRoot ? scopeRoot.querySelector('.calc-xyz-tray-controls') : null;
                        var trayStyle = tray ? window.getComputedStyle(tray) : null;
                        var trayVisible = !!(tray && trayStyle && trayStyle.display !== 'none');
                        var trayWidth = trayVisible ? tray.getBoundingClientRect().width : 0;
                        var availW = Math.max(120, rowRect.width - trayWidth - 16);
                        var h = Math.floor(availH * {CALC_MOL_DYNAMIC_SCALE});
                        var w = Math.floor(Math.min(h * 1.2, availW));
                        if (h >= 80 && w >= 120) {{
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
                {VIEWER_MOUSE_PATCH_JS}
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
                {VIEWER_MOUSE_PATCH_JS}
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
        large_traj_note = (
            f' <span style="color:#888;font-size:0.85em;">'
            f'(large trajectory, single-frame mode)</span>'
            if len(frames) > CALC_XYZ_LARGE_TRAJ_FRAMES else ''
        )
        calc_xyz_frame_total.value = f"<b>/ {len(frames)}</b>"
        calc_xyz_frame_label.value = (
            f"{format_xyz_comment_label(comment, frames=frames)}"
            f"{large_traj_note}"
        )
        selected_path = _calc_selected_item_path()
        rmsd_enabled = bool(
            selected_path
            and selected_path.suffix.lower() == '.xyz'
            and len(frames) == 1
        )
        _calc_set_rmsd_available(rmsd_enabled)
        _calc_update_xyz_traj_control_state()

        # Trajectory: use JS viewer and keep orientation when switching frames.
        # For large trajectories, fall back to single-frame mode to avoid
        # embedding tens of MB in a JS template literal over the comm channel.
        large_traj = len(frames) > CALC_XYZ_LARGE_TRAJ_FRAMES
        if len(frames) > 1 and not large_traj:
            if initial_load or not state['traj_viewer_ready']:
                full_xyz = ""
                for comm, block, natoms in frames:
                    full_xyz += f"{natoms}\n{comm}\n{block}\n"
                _mol3d_counter[0] += 1
                viewer_id = f"calc_trj_viewer_{_mol3d_counter[0]}"
                wrapper_id = f"calc_mol_wrap_{_mol3d_counter[0]}"
                traj_style_js = DEFAULT_3DMOL_STYLE_JS
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
                            /* Pre-size to actual free space in right panel */
                            var rightPanel = scopeRoot ? scopeRoot.querySelector('.calc-right') : null;
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
                                            if (!style || style.display === 'none' || style.visibility === 'hidden') {{
                                                continue;
                                            }}
                                            var childRect = child.getBoundingClientRect();
                                            if (childRect.height > 0) reservedBelow += childRect.height;
                                        }}
                                    }}
                                    var availH = rightRect.bottom - mvRect.top - reservedBelow - 10;
                                    var row = mv.closest('.calc-mol-view-row');
                                    var rowRect = row ? row.getBoundingClientRect() : rightRect;
                                    var tray = scopeRoot ? scopeRoot.querySelector('.calc-xyz-tray-controls') : null;
                                    var trayStyle = tray ? window.getComputedStyle(tray) : null;
                                    var trayVisible = !!(tray && trayStyle && trayStyle.display !== 'none');
                                    var trayWidth = trayVisible ? tray.getBoundingClientRect().width : 0;
                                    var availW = Math.max(120, rowRect.width - trayWidth - 16);
                                    var h = Math.floor(availH * {CALC_MOL_DYNAMIC_SCALE});
                                    var w = Math.floor(Math.min(h * 1.2, availW));
                                    if (h >= 80 && w >= 120) {{
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
                            {VIEWER_MOUSE_PATCH_JS}
                            var xyz = `{full_xyz}`;
                            viewer.addModelsAsFrames(xyz, "xyz");
                            viewer.setStyle({{}}, {traj_style_js});
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
                _calc_set_png_button_mode(main=True)
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
                # NOTE: _calc_update_png_frame was referenced here but never
                # defined — calling it would crash the PNG navigation button.
                # The JS above (v.setFrame(idx)) already updates the displayed
                # frame; the Python-side Python frame tracker isn't used.
                _calc_set_png_button_mode(main=True)
            return

        # Single frame: render with py3Dmol
        xyz_content = f"{n_atoms}\n{comment}\n{xyz_block}\n"
        _render_3dmol(xyz_content)
        _calc_apply_traj_style()

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

    def _calc_job_status_is_running(status):
        return str(status or '').upper().strip() in {
            'RUNNING', 'R', 'CG', 'COMPLETING',
        }

    def _calc_job_status_is_pending(status):
        return str(status or '').upper().strip() in {
            'PENDING', 'PD', 'Q', 'QUEUED', 'CF', 'CONFIGURING',
        }

    def _calc_job_status_activity_rank(status):
        if _calc_job_status_is_running(status):
            return 2
        if _calc_job_status_is_pending(status):
            return 1
        return 0

    def calc_collect_slurm_status_dirs(calc_root_resolved):
        """Collect SLURM statuses using workdir (%Z) for robust folder mapping."""
        if not shutil.which('squeue'):
            return []
        user = os.environ.get('USER')
        if not user:
            try:
                user = os.getlogin()
            except Exception:
                return []

        try:
            result = subprocess.run(
                ['squeue', '-h', '-u', user, '-o', '%i|%t|%Z'],
                capture_output=True,
                text=True,
                check=False,
                timeout=10,
            )
        except Exception:
            return []
        if result.returncode != 0:
            return []

        items = []
        for line in (result.stdout or '').splitlines():
            line = line.strip()
            if not line:
                continue
            parts = line.split('|', 2)
            if len(parts) < 3:
                continue
            raw_job_id, raw_status, raw_dir = parts
            raw_status = str(raw_status or '').upper().strip()
            raw_dir = str(raw_dir or '').strip()
            if not raw_status or not raw_dir:
                continue
            if raw_dir in {'N/A', '(null)', '-', 'None'}:
                continue

            try:
                job_id = int(str(raw_job_id).strip())
            except Exception:
                job_id = 0

            try:
                job_dir = Path(raw_dir).resolve()
                rel = job_dir.relative_to(calc_root_resolved)
            except Exception:
                continue
            if not rel.parts:
                continue
            folder_path = calc_root_resolved / rel.parts[0]
            items.append((folder_path, job_id, raw_status))
        return items

    def calc_collect_running_process_dirs():
        """Best-effort fallback: detect running DELFIN/ORCA process dirs."""
        user = os.environ.get('USER')
        if not user:
            try:
                user = os.getlogin()
            except Exception:
                return set()

        try:
            calc_root_resolved = _calc_dir().resolve()
        except Exception:
            calc_root_resolved = _calc_dir()

        keywords = ('orca', 'delfin', 'xtb', 'crest')
        try:
            result = subprocess.run(
                ['ps', '-u', user, '-o', 'pid=,command='],
                capture_output=True,
                text=True,
                check=False,
                timeout=10,
            )
        except Exception:
            return set()
        if result.returncode != 0:
            return set()

        running_root_dirs = set()
        for line in (result.stdout or '').splitlines():
            line = line.strip()
            if not line:
                continue
            parts = line.split(None, 1)
            if len(parts) < 2:
                continue
            pid, command = parts
            if not any(keyword in command.lower() for keyword in keywords):
                continue

            job_dir = ''
            try:
                job_dir = os.readlink(f'/proc/{pid}/cwd')
            except Exception:
                job_dir = ''

            if not job_dir:
                match = re.search(r'(/[^\\s]+\\.(?:inp|cisinp\\.tmp|goat|xyz))', command)
                if match:
                    try:
                        job_dir = str(Path(match.group(1)).parent)
                    except Exception:
                        job_dir = ''

            if not job_dir:
                match = re.search(r'(/[^\\s]+)', command)
                if match:
                    try:
                        p = Path(match.group(1))
                        if p.exists():
                            job_dir = str(p if p.is_dir() else p.parent)
                    except Exception:
                        job_dir = ''

            if not job_dir:
                continue

            try:
                resolved = Path(job_dir).resolve()
            except Exception:
                continue
            try:
                rel = resolved.relative_to(calc_root_resolved)
            except Exception:
                continue
            if not rel.parts:
                continue
            running_root_dirs.add(calc_root_resolved / rel.parts[0])

        return running_root_dirs

    def calc_collect_local_job_status_dirs():
        """Return latest local queue status per top-level calc folder."""
        backend = getattr(ctx, 'backend', None)
        load_jobs = getattr(backend, '_load_jobs', None)
        list_jobs = getattr(backend, 'list_jobs', None)
        data = {}
        if callable(load_jobs):
            try:
                loaded = load_jobs()
                if isinstance(loaded, dict):
                    data = loaded
            except Exception:
                data = {}

        calc_root = _calc_dir()
        try:
            calc_root_resolved = calc_root.resolve()
        except Exception:
            calc_root_resolved = calc_root

        latest_by_folder = {}
        def _update_latest(folder_path, job_id, raw_status):
            if folder_path is None:
                return
            prev = latest_by_folder.get(folder_path)
            if prev is None:
                latest_by_folder[folder_path] = (job_id, raw_status)
                return
            prev_id, prev_status = prev
            prev_rank = _calc_job_status_activity_rank(prev_status)
            now_rank = _calc_job_status_activity_rank(raw_status)
            if now_rank > prev_rank:
                latest_by_folder[folder_path] = (max(prev_id, job_id), raw_status)
            elif now_rank == prev_rank and job_id >= prev_id:
                latest_by_folder[folder_path] = (job_id, raw_status)

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

            _update_latest(folder_path, job_id, raw_status)

        # SLURM source with reliable workdir mapping (incl. PD pending jobs).
        for folder_path, job_id, raw_status in calc_collect_slurm_status_dirs(calc_root_resolved):
            _update_latest(folder_path, job_id, raw_status)

        # Also consider live backend jobs (e.g., SLURM list_jobs output),
        # using job_dir when available and falling back to calc/<job_name>.
        if callable(list_jobs):
            try:
                live_jobs = list_jobs() or []
            except Exception:
                live_jobs = []
            for job in live_jobs:
                raw_status = str(getattr(job, 'status', '')).upper().strip()
                if not raw_status:
                    continue

                folder_path = None
                raw_dir = str(getattr(job, 'job_dir', '') or '').strip()
                if raw_dir:
                    try:
                        job_dir = Path(raw_dir).resolve()
                        rel = job_dir.relative_to(calc_root_resolved)
                        if rel.parts:
                            folder_path = calc_root_resolved / rel.parts[0]
                    except Exception:
                        folder_path = None

                if folder_path is None:
                    job_name = str(getattr(job, 'name', '') or '').strip()
                    if job_name:
                        try:
                            candidate = (calc_root_resolved / job_name).resolve()
                        except Exception:
                            candidate = calc_root_resolved / job_name
                        try:
                            candidate.relative_to(calc_root_resolved)
                        except Exception:
                            candidate = None
                        if candidate is not None and candidate.exists() and candidate.is_dir():
                            folder_path = candidate

                if folder_path is None:
                    continue

                raw_job_id = getattr(job, 'job_id', 0)
                try:
                    job_id = int(str(raw_job_id))
                except Exception:
                    job_id = 0

                _update_latest(folder_path, job_id, raw_status)

        # Overlay live process detection so running jobs still get a blue dot
        # when the queue metadata is temporarily empty/stale.
        for folder_path in calc_collect_running_process_dirs():
            prev = latest_by_folder.get(folder_path)
            if prev is None:
                latest_by_folder[folder_path] = (0, 'RUNNING')
            else:
                latest_by_folder[folder_path] = (prev[0], 'RUNNING')

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
        _calc_stop_xyz_playback(update_button=True)
        state['xyz_frames'].clear()
        state['xyz_current_frame'][0] = 0
        _calc_clear_main_viewer_state(reset_view_state=False)
        calc_mol_container.layout.display = 'none'
        calc_xyz_controls.layout.display = 'none'
        calc_coord_controls.layout.display = 'none'
        _calc_update_xyz_traj_control_state()
        _calc_set_rmsd_available(False)
        _calc_reset_rmsd_ui(clear_input=True)
        calc_file_info.value = ''
        calc_path_display.value = ''
        calc_download_status.value = ''
        calc_xyz_batch_status.value = ''
        state['xyz_batch_marked_paths'] = []
        calc_print_mode_status.value = ''
        calc_mo_plot_status.value = ''
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
        _calc_show_xyz_batch_panel(False)
        _calc_show_print_mode_panel(False)
        _calc_show_mo_plot_panel(False)
        _calc_hide_chunk_controls()
        calc_update_view()
        calc_set_message('Select a file...')
        _calc_process_staged_uploads()

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
                is_top_level_calc_view = (
                    not _is_archive_tab
                    and current_dir.resolve() == _calc_dir().resolve()
                )
            except Exception:
                is_top_level_calc_view = (
                    not _is_archive_tab
                    and (state.get('current_path') or '') == ''
                )
            local_status_dirs = calc_collect_local_job_status_dirs()
            for entry in entries:
                if entry.is_dir():
                    if is_top_level_calc_view:
                        is_running = False
                        is_pending = False
                        try:
                            entry_resolved = entry.resolve()
                        except Exception:
                            entry_resolved = entry
                        if entry_resolved in local_status_dirs:
                            status = local_status_dirs.get(entry_resolved, '')
                            is_running = _calc_job_status_is_running(status)
                            is_pending = _calc_job_status_is_pending(status)
                            folder_icon = '✅' if status == 'COMPLETED' else '📂'
                        elif calc_folder_has_completed_results(entry):
                            folder_icon = '✅'
                        else:
                            folder_icon = '📂'
                        if is_running:
                            folder_icon = '🔵'
                        elif is_pending:
                            folder_icon = '🟠'
                    else:
                        folder_icon = '📂'
                    items.append(f'{folder_icon} {entry.name}')
                else:
                    suffix = entry.suffix.lower()
                    if _calc_is_table_preset_path(entry):
                        items.append(f'📊 {entry.name}')
                    elif suffix == '.xyz':
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
        calc_update_download_btn()
        _calc_update_explorer_action_state()

    def calc_filter_file_list(change=None):
        query = calc_folder_search.value.strip().lower()
        if not query:
            calc_file_list.options = state['all_items']
        else:
            filtered = [item for item in state['all_items'] if query in item.lower()]
            result = filtered if filtered else ['(No matches)']
            calc_file_list.options = result

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

    def calc_reset_xyz_workflow_state():
        state['xyz_workflow_active'] = False
        state['xyz_workflow_mode'] = ''
        state['xyz_workflow_source_path'] = None
        calc_xyz_workflow_pal.value = CALC_BROWSER_WORKFLOW_PAL
        calc_xyz_workflow_time.value = CALC_BROWSER_WORKFLOW_TIMELIMIT
        calc_xyz_workflow_bfw.value = False
        calc_xyz_workflow_bfw.description = 'BFW Off'
        calc_xyz_workflow_bfw.button_style = ''
        calc_xyz_workflow_bfw.layout.display = 'none'
        calc_submit_xyz_workflow_btn.disabled = True
        calc_xyz_workflow_info.value = ''
        calc_xyz_workflow_status.value = ''
        calc_xyz_workflow_toolbar.layout.display = 'none'

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
        candidates = []
        for i in range(1, len(path_parts) + 1):
            test_path = '/'.join(path_parts[:i])
            test_dir = _calc_dir() / test_path
            if (test_dir / 'CONTROL.txt').exists():
                candidates.append(test_path)

        if not candidates:
            return None

        # Prefer the outermost real DELFIN workspace (contains top-level markers),
        # not nested OCCUPIER subfolders that also carry CONTROL.txt snapshots.
        for candidate in candidates:
            candidate_dir = _calc_dir() / candidate
            if (
                (candidate_dir / 'input.txt').exists()
                or (candidate_dir / 'start.txt').exists()
                or (candidate_dir / 'DELFIN.txt').exists()
            ):
                return candidate

        # Fallback: still use the outermost CONTROL-bearing directory.
        return candidates[0]

    def _calc_is_report_workspace(path):
        if not path.is_dir() or not (path / 'CONTROL.txt').exists():
            return False
        return any(
            (path / name).exists()
            for name in (
                'input.txt',
                'start.txt',
                'DELFIN.txt',
                'DELFIN_Data.json',
                'DELFIN.docx',
                'report.docx',
            )
        )

    def _calc_report_target_label(path):
        try:
            return str(path.relative_to(_calc_dir()))
        except Exception:
            return path.name or str(path)

    def calc_collect_report_targets(limit=None):
        current_dir = _calc_dir() / state['current_path'] if state['current_path'] else _calc_dir()
        workspace_root = calc_find_workspace_root(state['current_path']) if state['current_path'] else None
        if workspace_root:
            workspace_dir = _calc_dir() / workspace_root
            if _calc_is_report_workspace(workspace_dir):
                return [workspace_dir]

        targets = []
        seen = set()

        def _add_target(path):
            try:
                resolved = path.resolve()
            except Exception:
                resolved = path
            if resolved in seen:
                return False
            seen.add(resolved)
            targets.append(path)
            return True

        if _calc_is_report_workspace(current_dir):
            _add_target(current_dir)
            if limit is not None and len(targets) >= limit:
                return targets

        try:
            for control_path in current_dir.rglob('CONTROL.txt'):
                candidate = control_path.parent
                if not _calc_is_report_workspace(candidate):
                    continue
                if _add_target(candidate) and limit is not None and len(targets) >= limit:
                    break
        except Exception:
            pass

        try:
            targets.sort(
                key=lambda path: (
                    len(path.relative_to(current_dir).parts),
                    str(path).lower(),
                )
            )
        except Exception:
            targets.sort(key=lambda path: str(path).lower())

        return targets

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
        has_multiple = len(state['search_spans']) > 1
        calc_prev_btn.disabled = not has_multiple
        calc_next_btn.disabled = not has_multiple

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
        else:
            start, end = state['search_spans'][state['current_match']]
            _calc_focus_plain_match(start, end)

    # -- options dropdown ---------------------------------------------------
    def calc_update_options_dropdown():
        labels = _calc_selected_labels()
        selected = labels[0] if labels else ''
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
        if selected and sel_lower.endswith('.xyz'):
            xyz_options = ['(Options)', 'Build Batch from XYZ', 'Calc NMR', 'Calc CENSO/ANMR']
            selected_path = _calc_selected_item_path()
            if selected_path and _calc_is_single_structure_xyz(selected_path):
                xyz_options.extend(['hyperpol_xtb', 'tadf_xtb'])
            if rmsd_available:
                xyz_options.append('RMSD')
            calc_options_dropdown.options = xyz_options
            calc_options_dropdown.value = '(Options)'
            calc_options_dropdown.layout.display = 'block'
            return
        if selected and 'OCCUPIER.txt' in selected:
            calc_options_dropdown.options = ['(Options)', 'Override']
            calc_options_dropdown.value = '(Options)'
            calc_options_dropdown.layout.display = 'block'
        elif selected and 'CONTROL.txt' in selected:
            calc_options_dropdown.options = ['(Options)', 'Recalc', 'Smart Recalc']
            calc_options_dropdown.value = '(Options)'
            calc_options_dropdown.layout.display = 'block'
        elif selected and sel_lower.endswith('.inp'):
            calc_options_dropdown.options = ['(Options)', 'Recalc']
            calc_options_dropdown.value = '(Options)'
            calc_options_dropdown.layout.display = 'block'
        elif selected and sel_lower.endswith('.out'):
            out_options = ['(Options)', 'Print Mode', 'MO Plot', 'Print NMR']
            if rmsd_available:
                out_options.append('RMSD')
            calc_options_dropdown.options = out_options
            calc_options_dropdown.value = '(Options)'
            calc_options_dropdown.layout.display = 'block'
        elif selected and sel_lower.endswith('.final.interp'):
            calc_options_dropdown.options = ['(Options)', 'Plot Trajectory']
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
            if (
                not calc_view_toggle.value
                and not state['recalc_active']
                and not state['xyz_workflow_active']
            ):
                calc_content_area.layout.display = 'block'

    # -- event handlers -----------------------------------------------------
    def calc_on_back(b):
        if state['current_path']:
            _calc_hide_rename_prompt()
            parts = state['current_path'].split('/')
            state['current_path'] = '/'.join(parts[:-1])
            calc_list_directory()
            calc_file_list.value = ()

    def calc_on_home(b):
        _calc_hide_rename_prompt()
        state['current_path'] = ''
        calc_list_directory()
        calc_file_list.value = ()

    def calc_on_refresh(b):
        _calc_hide_rename_prompt()
        calc_list_directory()
        calc_file_list.value = ()

    def _calc_resolve_navigation_target(raw_value):
        raw = str(raw_value or '').strip()
        root = _calc_dir().resolve()
        if raw in {'', '/'}:
            return root
        if raw.startswith(str(root)):
            candidate = Path(raw).expanduser()
        elif raw.startswith('/'):
            candidate = root / raw.lstrip('/')
        else:
            candidate = _calc_current_dir() / raw
        return _calc_resolve_within_root(candidate)

    def calc_on_path_change(change):
        if _calc_path_syncing[0]:
            return
        raw = str(change.get('new') or '').strip()
        if not raw:
            calc_update_path_display()
            return
        try:
            target = _calc_resolve_navigation_target(raw)
            if not target.exists():
                raise FileNotFoundError('Path not found inside explorer root.')
            calc_folder_search.value = ''
            if target.is_dir():
                state['current_path'] = _calc_state_rel_path(target)
                calc_list_directory()
                calc_file_list.value = ()
            else:
                state['current_path'] = _calc_state_rel_path(target.parent)
                calc_list_directory()
                target_label = _calc_label_for_name(target.name)
                if not target_label:
                    raise FileNotFoundError('Target file could not be selected in explorer.')
                calc_file_list.value = (target_label,)
                _calc_open_item(target_label)
            _calc_set_ops_status('', '#555')
        except Exception as exc:
            calc_update_path_display()
            _calc_set_ops_status(_html.escape(str(exc)), '#d32f2f')

    def calc_on_back_to_calculations(_button=None):
        # Archive-tab shortcut: same behavior as move-button, but opposite direction.
        calc_on_move_archive_click(_button)

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
            else:
                _calc_clear_plain_match_marker()
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
            else:
                _calc_clear_plain_match_marker()
            return

        state['current_match'] = 0
        calc_update_nav_buttons()
        calc_update_search_result()
        if _calc_should_highlight():
            calc_apply_highlight(query, state['current_match'])
        else:
            start, end = state['search_spans'][state['current_match']]
            _calc_focus_plain_match(start, end)

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
        if not state['search_spans']:
            return
        if len(state['search_spans']) == 1:
            state['current_match'] = 0
            calc_show_match()
            return
        state['current_match'] = (
            state['current_match'] - 1
            if state['current_match'] > 0
            else len(state['search_spans']) - 1
        )
        calc_show_match()

    def calc_next_match(b):
        if not state['search_spans']:
            return
        if len(state['search_spans']) == 1:
            state['current_match'] = 0
            calc_show_match()
            return
        state['current_match'] = (
            state['current_match'] + 1
            if state['current_match'] < len(state['search_spans']) - 1
            else 0
        )
        calc_show_match()

    def calc_on_xyz_input_change(change):
        new_val = change['new']
        if state['xyz_frames'] and 1 <= new_val <= len(state['xyz_frames']):
            state['xyz_current_frame'][0] = new_val - 1
            calc_update_xyz_viewer(initial_load=False)

    def calc_on_xyz_loop_change(change):
        if change.get('name') != 'value':
            return
        _calc_update_loop_button_style()

    def calc_on_xyz_fps_change(change):
        if change.get('name') != 'value':
            return
        try:
            fps = int(change.get('new'))
        except Exception:
            fps = CALC_XYZ_PLAY_FPS_DEFAULT
        fps = max(CALC_XYZ_PLAY_FPS_MIN, min(CALC_XYZ_PLAY_FPS_MAX, fps))
        if calc_xyz_fps_input.value != fps:
            calc_xyz_fps_input.value = fps

    def calc_on_xyz_play_change(change):
        if change.get('name') != 'value':
            return
        if state.get('traj_play_toggle_guard'):
            return
        if bool(change.get('new')):
            _calc_start_xyz_playback()
        else:
            _calc_stop_xyz_playback(update_button=True)

    def calc_on_xyz_copy(button):
        frames = state['xyz_frames']
        if not frames:
            return
        idx = state['xyz_current_frame'][0]
        if idx < 0 or idx >= len(frames):
            return
        comment, xyz_block, n_atoms = frames[idx]
        try:
            symbols, coords = _calc_parse_xyz_lines(
                xyz_block.splitlines(), expected_atoms=n_atoms, allow_skip=False
            )
            xyz_content = _calc_build_xyz_from_symbols_coords(
                symbols, coords, comment=(comment or f'Frame {idx + 1}')
            ).rstrip()
        except Exception:
            xyz_content = f"{n_atoms}\n{comment}\n{xyz_block}"
        _calc_copy_to_clipboard(xyz_content, label=f'xyz frame {idx + 1}')

    def calc_on_view_png(button):
        _calc_trigger_png_download()

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
        _calc_copy_to_clipboard(state['file_content'], label='coord content')

    def calc_on_content_copy(button):
        content = ''
        if state['recalc_active']:
            content = calc_edit_area.value or ''
        if not content:
            content = state['file_content'] or ''
        if not content:
            return
        _calc_copy_to_clipboard(content, label='file content')

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
        _calc_copy_to_clipboard(full_path, label='path')

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

        # If folders are selected, only generate reports for those.
        selected = _calc_selected_labels()
        selected_dirs = []
        if selected:
            for label in selected:
                name = _calc_label_to_name(label)
                if not name:
                    continue
                candidate = current_dir / name
                if candidate.is_dir() and _calc_is_report_workspace(candidate):
                    selected_dirs.append(candidate)
                elif candidate.is_dir():
                    # Check subfolders of selected dir
                    try:
                        for ctrl in candidate.rglob('CONTROL.txt'):
                            sub = ctrl.parent
                            if _calc_is_report_workspace(sub):
                                selected_dirs.append(sub)
                    except Exception:
                        pass

        report_targets = selected_dirs if selected_dirs else calc_collect_report_targets()
        if not report_targets:
            calc_report_status.value = (
                '<span style="color:#d32f2f;">No reportable DELFIN folders found here.</span>'
            )
            return

        dir_key = str(current_dir)

        if dir_key in state['report_running']:
            calc_report_status.value = '<span style="color:orange;">Report already running for this folder...</span>'
            return

        state['report_running'][dir_key] = True
        if len(report_targets) == 1:
            target_label = _calc_report_target_label(report_targets[0])
            calc_report_status.value = (
                f'<span style="color:blue;">Generating report in {target_label}...</span>'
            )
        else:
            calc_report_status.value = (
                f'<span style="color:blue;">Generating reports in {len(report_targets)} folders...</span>'
            )

        def run_report():
            successes = []
            failures = []
            timeouts = []
            try:
                for idx, target_dir in enumerate(report_targets, start=1):
                    target_label = _calc_report_target_label(target_dir)
                    if len(report_targets) == 1:
                        calc_report_status.value = (
                            f'<span style="color:blue;">Generating report in {target_label}...</span>'
                        )
                    else:
                        calc_report_status.value = (
                            f'<span style="color:blue;">Generating report {idx}/{len(report_targets)} '
                            f'in {target_label}...</span>'
                        )

                    try:
                        result = subprocess.run(
                            ['delfin', '--report', 'docx'],
                            cwd=str(target_dir),
                            capture_output=True, text=True, timeout=600,
                        )
                        if result.returncode == 0:
                            successes.append(target_label)
                        else:
                            msg = (result.stderr or result.stdout or 'Unknown error').strip()
                            failures.append((target_label, msg[:160]))
                    except subprocess.TimeoutExpired:
                        timeouts.append(target_label)
                    except Exception as exc:
                        failures.append((target_label, str(exc)[:160]))

                if failures or timeouts:
                    parts = []
                    if successes:
                        parts.append(f'{len(successes)} ok')
                    if failures:
                        parts.append(f'{len(failures)} failed')
                    if timeouts:
                        parts.append(f'{len(timeouts)} timeout')
                    detail = ''
                    if failures:
                        detail = f' First error in {failures[0][0]}: {_html.escape(failures[0][1])}'
                    elif timeouts:
                        detail = f' Timeout in {timeouts[0]}.'
                    calc_report_status.value = (
                        f'<span style="color:#d32f2f;">Report run finished: {", ".join(parts)}.'
                        f'{detail}</span>'
                    )
                else:
                    if len(successes) == 1:
                        calc_report_status.value = (
                            f'<span style="color:green;">Report generated in {successes[0]}!</span>'
                        )
                    else:
                        calc_report_status.value = (
                            f'<span style="color:green;">Reports generated in {len(successes)} folders!</span>'
                        )
            except Exception as e:
                calc_report_status.value = (
                    f'<span style="color:red;">Error: {_html.escape(str(e)[:100])}</span>'
                )
            finally:
                state['report_running'].pop(dir_key, None)

        threading.Thread(target=run_report, daemon=True).start()

    def calc_on_options_change(change):
        if change['new'] != 'Build Batch from XYZ':
            _calc_show_xyz_batch_panel(False)
        if change['new'] != 'Print Mode':
            _calc_show_print_mode_panel(False)
        if change['new'] != 'MO Plot':
            _calc_show_mo_plot_panel(False)
        if change['new'] != 'Calc NMR':
            _calc_show_nmr_panel(False)
        if change['new'] != 'Calc CENSO/ANMR':
            _calc_show_censo_nmr_panel(False)
        if change['new'] not in ('hyperpol_xtb', 'tadf_xtb'):
            calc_reset_xyz_workflow_state()
        if change['new'] == 'Override':
            calc_override_input.layout.display = 'block'
            calc_override_time.layout.display = 'block'
            calc_override_btn.layout.display = 'block'
            calc_override_status.layout.display = 'block'
            calc_override_status.value = ''
            calc_edit_area.layout.display = 'none'
            _calc_preselect_show(False)
            calc_content_area.layout.display = 'block'
            cur_labels = _calc_selected_labels()
            if cur_labels:
                name = _calc_label_to_name(cur_labels[0])
                stage_name = name.replace('_OCCUPIER.txt', '').replace('OCCUPIER.txt', '')
                calc_override_input.value = f'{stage_name}=' if stage_name else ''
        elif change['new'] in ('Recalc', 'Smart Recalc'):
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
        elif change['new'] == 'Build Batch from XYZ':
            calc_override_input.layout.display = 'none'
            calc_override_time.layout.display = 'none'
            calc_override_btn.layout.display = 'none'
            calc_override_status.layout.display = 'none'
            calc_override_status.value = ''
            calc_edit_area.layout.display = 'none'
            _calc_preselect_show(False)
            selected_path = _calc_selected_item_path()
            if not selected_path or selected_path.suffix.lower() != '.xyz':
                _calc_show_xyz_batch_panel(False)
                calc_xyz_batch_status.value = (
                    '<span style="color:#d32f2f;">Select a .xyz file first.</span>'
                )
                return
            state['xyz_batch_marked_paths'] = []
            state['xyz_batch_path'] = ''
            calc_xyz_batch_filename.value = _calc_xyz_batch_default_filename()
            _calc_show_xyz_batch_panel(True)
        elif change['new'] in ('hyperpol_xtb', 'tadf_xtb'):
            calc_override_input.layout.display = 'none'
            calc_override_time.layout.display = 'none'
            calc_override_btn.layout.display = 'none'
            calc_override_status.layout.display = 'none'
            calc_override_status.value = ''
            calc_edit_area.layout.display = 'none'
            _calc_preselect_show(False)
            _calc_show_xyz_batch_panel(False)
            calc_content_area.layout.display = 'block'
            _calc_prepare_xyz_browser_workflow(change['new'])
            _calc_reset_options_dropdown()
        elif change['new'] == 'Print Mode':
            calc_override_input.layout.display = 'none'
            calc_override_time.layout.display = 'none'
            calc_override_btn.layout.display = 'none'
            calc_override_status.layout.display = 'none'
            calc_override_status.value = ''
            calc_edit_area.layout.display = 'none'
            _calc_preselect_show(False)
            _calc_show_xyz_batch_panel(False)
            _calc_show_print_mode_panel(True)
            calc_content_area.layout.display = 'block'
        elif change['new'] == 'MO Plot':
            calc_override_input.layout.display = 'none'
            calc_override_time.layout.display = 'none'
            calc_override_btn.layout.display = 'none'
            calc_override_status.layout.display = 'none'
            calc_override_status.value = ''
            calc_edit_area.layout.display = 'none'
            _calc_preselect_show(False)
            _calc_show_xyz_batch_panel(False)
            _calc_show_print_mode_panel(False)
            _calc_show_mo_plot_panel(True)
            calc_content_area.layout.display = 'block'
        elif change['new'] == 'Plot Trajectory':
            calc_override_input.layout.display = 'none'
            calc_override_time.layout.display = 'none'
            calc_override_btn.layout.display = 'none'
            calc_override_status.layout.display = 'none'
            calc_override_status.value = ''
            calc_edit_area.layout.display = 'none'
            _calc_preselect_show(False)
            _calc_show_xyz_batch_panel(False)
            _calc_show_print_mode_panel(False)
            _calc_show_mo_plot_panel(False)
            _calc_render_traj_plot()
            calc_content_area.layout.display = 'block'
            _calc_reset_options_dropdown()
        elif change['new'] == 'Print NMR':
            calc_override_input.layout.display = 'none'
            calc_override_time.layout.display = 'none'
            calc_override_btn.layout.display = 'none'
            calc_override_status.layout.display = 'none'
            calc_override_status.value = ''
            calc_edit_area.layout.display = 'none'
            _calc_preselect_show(False)
            _calc_show_xyz_batch_panel(False)
            _calc_show_print_mode_panel(False)
            _calc_show_mo_plot_panel(False)
            calc_content_area.layout.display = 'block'
            _calc_run_print_nmr()
            _calc_reset_options_dropdown()
        elif change['new'] == 'Calc NMR':
            calc_override_input.layout.display = 'none'
            calc_override_time.layout.display = 'none'
            calc_override_btn.layout.display = 'none'
            calc_override_status.layout.display = 'none'
            calc_override_status.value = ''
            calc_edit_area.layout.display = 'none'
            _calc_preselect_show(False)
            _calc_show_xyz_batch_panel(False)
            _calc_show_print_mode_panel(False)
            _calc_show_mo_plot_panel(False)
            _calc_show_nmr_panel(True)
            calc_content_area.layout.display = 'block'
        elif change['new'] == 'Calc CENSO/ANMR':
            calc_override_input.layout.display = 'none'
            calc_override_time.layout.display = 'none'
            calc_override_btn.layout.display = 'none'
            calc_override_status.layout.display = 'none'
            calc_override_status.value = ''
            calc_edit_area.layout.display = 'none'
            _calc_preselect_show(False)
            _calc_show_xyz_batch_panel(False)
            _calc_show_print_mode_panel(False)
            _calc_show_mo_plot_panel(False)
            _calc_show_censo_nmr_panel(True)
            calc_content_area.layout.display = 'block'
        elif change['new'] == 'RMSD':
            calc_override_input.layout.display = 'none'
            calc_override_time.layout.display = 'none'
            calc_override_btn.layout.display = 'none'
            calc_override_status.layout.display = 'none'
            calc_override_status.value = ''
            calc_edit_area.layout.display = 'none'
            _calc_preselect_show(False)
            _calc_show_xyz_batch_panel(False)
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
            _calc_show_xyz_batch_panel(False)
            calc_content_area.layout.display = 'block'

    def calc_on_override_start(button):
        selected_option = calc_options_dropdown.value
        is_recalc = selected_option in ('Recalc', 'Smart Recalc')
        is_smart_recalc = selected_option == 'Smart Recalc'
        _ovr_labels = _calc_selected_labels()
        selected_name = _calc_label_to_name(_ovr_labels[0]) if _ovr_labels else ''
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
                recalc_mode = 'delfin-recalc' if is_smart_recalc else 'delfin-recalc-classic'
                result = ctx.backend.submit_delfin(
                    job_dir=workspace_path,
                    job_name=job_name,
                    mode=recalc_mode,
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
                if not is_recalc:
                    calc_override_status.value = (
                        f'<span style="color:#388e3c;">Submitted from {workspace_root} '
                        f'(Override: {_html.escape(override_value)}): '
                        f'{result.stdout.strip()}</span>'
                    )
                else:
                    recalc_label = 'Smart Recalc' if is_smart_recalc else 'Recalc'
                    calc_override_status.value = (
                        f'<span style="color:#388e3c;">Submitted {recalc_label} from {workspace_root}: '
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

    def calc_on_delete_click(button=None):
        real = _calc_selected_labels()
        if real:
            names = [_calc_label_to_name(s) for s in real]
            label_str = ', '.join(f'<code>{_html.escape(n)}</code>' for n in names)
            calc_delete_label.value = f'<b>Delete {len(names)} item(s)?</b> {label_str}'
            state['delete_current'] = False
            calc_delete_confirm.layout.display = 'flex'
            return
        if state['current_path']:
            calc_delete_label.value = (
                f'<b>Delete current folder?</b> '
                f'<code>{_html.escape(state["current_path"])}</code>'
            )
            state['delete_current'] = True
            calc_delete_confirm.layout.display = 'flex'
        else:
            calc_delete_status.value = '<span style="color:#d32f2f;">No item selected.</span>'

    def calc_on_delete_yes(button):
        if not state.get('delete_current', False):
            real = _calc_selected_labels()
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
            calc_file_list.value = ()
            _calc_refresh_related_explorers()
            return
        # delete_current fallback (delete the folder we are inside)
        if not state['current_path']:
            calc_delete_status.value = '<span style="color:#d32f2f;">No current folder.</span>'
            calc_delete_confirm.layout.display = 'none'
            state['delete_current'] = False
            return
        full_path = _calc_dir() / state['current_path']
        name = state['current_path']
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
        # This code path IS the "delete current folder" fallback (see comment
        # above at line 9191), so after successful deletion we always move up
        # to the parent directory. Previously this guarded on an undefined
        # local `delete_current` (a leftover after a state-dict refactor);
        # the state-dict flag was cleared before the guard ran, so the guard
        # body never executed and the parent-navigation was silently dead.
        state['delete_current'] = False
        cp = state['current_path']
        state['current_path'] = cp.rsplit('/', 1)[0] if '/' in cp else ''
        calc_file_list.value = ()
        _calc_refresh_related_explorers()

    def calc_on_delete_no(button):
        calc_delete_hide_confirm()

    def _calc_collect_archive_toggle_sources():
        labels = _calc_selected_labels()
        current_dir = _calc_current_dir()
        sources = []
        for label in labels:
            name = _calc_label_to_name(label)
            if not name:
                continue
            src = current_dir / name
            if src.exists():
                sources.append(src)
        if sources:
            return sources
        # Fallback requested by users: when nothing is selected, move current folder.
        if state['current_path']:
            cur = _calc_current_dir()
            if cur.exists() and cur.resolve() != _calc_dir().resolve():
                return [cur]
        return []

    def _calc_archive_toggle_destination():
        if _is_archive_tab:
            calc_root = getattr(ctx, 'primary_calc_dir', None)
            if not calc_root:
                calc_root = ctx.archive_dir.parent / 'calc'
            return Path(calc_root), 'Calculations'
        return ctx.archive_dir, 'Archive'

    def calc_on_move_archive_click(button):
        sources = _calc_collect_archive_toggle_sources()
        if not sources:
            calc_move_archive_status.value = (
                '<span style="color:#d32f2f;">No items selected and no current folder to move.</span>'
            )
            return
        dest_dir, dest_label = _calc_archive_toggle_destination()
        names = [src.name for src in sources]
        label_str = ', '.join(f'<code>{_html.escape(n)}</code>' for n in names)
        calc_move_archive_label.value = (
            f'<b>Move {len(names)} item(s) to {dest_label}?</b> {label_str}'
        )
        calc_move_archive_confirm.layout.display = 'flex'

    def calc_on_move_archive_yes(button):
        sources = _calc_collect_archive_toggle_sources()
        calc_move_archive_confirm.layout.display = 'none'
        if not sources:
            calc_move_archive_status.value = (
                '<span style="color:#d32f2f;">No items selected and no current folder to move.</span>'
            )
            return

        dest_dir, dest_label = _calc_archive_toggle_destination()
        dest_dir = Path(dest_dir)
        dest_dir.mkdir(parents=True, exist_ok=True)
        dest_resolved = dest_dir.resolve()

        errors, moved = [], []
        current_before = _calc_current_dir().resolve()
        moved_current_folder = False

        for src in sources:
            name = src.name
            try:
                src_resolved = _calc_resolve_within_root(src)
                dst = dest_resolved / src_resolved.name
                if dst.exists():
                    raise ValueError(f'Target already exists: {src_resolved.name}')
                shutil.move(str(src_resolved), str(dst))
                moved.append(src_resolved.name)
                if src_resolved == current_before:
                    moved_current_folder = True
            except Exception as exc:
                errors.append(f'{_html.escape(name)}: {_html.escape(str(exc))}')

        if moved_current_folder:
            cp = state.get('current_path', '')
            state['current_path'] = cp.rsplit('/', 1)[0] if '/' in cp else ''

        if moved and not errors:
            calc_move_archive_status.value = (
                f'<span style="color:#2e7d32;">Moved {len(moved)} item(s) to {dest_label}.</span>'
            )
        elif moved and errors:
            calc_move_archive_status.value = (
                f'<span style="color:#e65100;">Moved {len(moved)} to {dest_label}, '
                f'failed: {"; ".join(errors)}</span>'
            )
        else:
            calc_move_archive_status.value = (
                f'<span style="color:#d32f2f;">{"; ".join(errors)}</span>'
            )

        calc_file_list.value = ()
        calc_file_list.value = ()
        _calc_refresh_related_explorers()

    def calc_on_move_archive_no(button):
        calc_move_archive_confirm.layout.display = 'none'

    def _calc_start_ssh_transfer_with_config(sources, config_payload):
        if state.get('ssh_transfer_running', False):
            return
        resolved_sources = [_calc_resolve_within_root(src) for src in sources]
        host, user, remote_path, port = normalize_ssh_transfer_settings(
            config_payload.get('host', ''),
            config_payload.get('user', ''),
            config_payload.get('remote_path', ''),
            config_payload.get('port', 22),
        )
        job = create_transfer_job(
            resolved_sources,
            host,
            user,
            remote_path,
            port,
            delete_local_on_success=True,
        )
        source_count = len(resolved_sources)
        target_label = _html.escape(f'{user}@{host}:{remote_path}')
        _calc_set_ssh_transfer_running(True)
        try:
            pid = launch_transfer_job(job, python_executable=sys.executable)
            calc_transfer_jobs_panel.layout.display = 'flex'
            _calc_update_transfer_jobs_visibility()
            _calc_render_transfer_jobs()
            _calc_set_ops_status(
                f'Started resumable transfer job <code>{_html.escape(job["job_id"])}</code> '
                f'for {source_count} item(s) to <code>{target_label}</code>. '
                'Local sources will be removed only after rsync finishes successfully. '
                f'PID: <code>{int(pid)}</code>. Log: <code>{_html.escape(job["log_path"])}</code>.',
                color='#2e7d32',
            )
        except Exception as exc:
            _calc_set_ops_status(
                f'SSH transfer failed: {_html.escape(str(exc))}',
                color='#d32f2f',
            )
        finally:
            _calc_set_ssh_transfer_running(False)

    def calc_on_ssh_transfer_click(button=None):
        _calc_hide_rename_prompt()
        sources = _calc_collect_selected_sources_only()
        if not sources:
            _calc_update_explorer_action_state()
            _calc_set_ops_status('Select item(s) first.', color='#d32f2f')
            return

        config_payload = _calc_load_transfer_config()
        if config_payload is False:
            return
        if not config_payload:
            ctx.select_tab('Settings')
            _calc_set_ops_status(
                (
                    f'No transfer target saved yet. Opened <b>Settings</b>. '
                    f'Save host, user, port and target once in '
                    f'<code>{_html.escape(str(_calc_settings_path()))}</code>, '
                    'then the SSH Transfer button can start directly.'
                ),
                color='#ef6c00',
            )
            return
        _calc_start_ssh_transfer_with_config(sources, config_payload)

    def calc_on_transfer_jobs_toggle(button=None):
        if calc_transfer_jobs_panel.layout.display == 'none':
            _calc_hide_rename_prompt()
            _calc_render_transfer_jobs()
            calc_transfer_jobs_panel.layout.display = 'flex'
        else:
            calc_transfer_jobs_panel.layout.display = 'none'
        _calc_update_transfer_jobs_visibility()

    def calc_on_transfer_jobs_refresh(button=None):
        _calc_render_transfer_jobs()

    def _calc_parse_move_target_dir(raw_target):
        target = str(raw_target or '').strip()
        if not target:
            raise ValueError('Provide destination folder path.')
        if target.startswith('/'):
            rel = target.lstrip('/')
            target_path = (_calc_dir() / rel) if rel else _calc_dir()
        else:
            target_path = _calc_current_dir() / target
        target_path = _calc_resolve_within_root(target_path)
        if not target_path.exists() or not target_path.is_dir():
            raise ValueError(f'Destination folder not found: {target}')
        return target_path

    def calc_on_new_folder(button=None):
        try:
            raw_name = str(calc_new_folder_input.value or '').strip()
            folder_name = (
                _calc_clean_item_name(raw_name)
                if raw_name
                else _calc_next_available_folder_name('New Folder')
            )
            target = _calc_resolve_within_root(_calc_current_dir() / folder_name)
            if target.exists():
                raise ValueError(f'Folder already exists: {folder_name}')
            target.mkdir(parents=False, exist_ok=False)
            calc_new_folder_input.value = ''
            _calc_set_ops_status(
                f'Folder created: <code>{_html.escape(folder_name)}</code>',
                color='#2e7d32',
            )
            calc_file_list.value = ()
            _calc_refresh_related_explorers()
        except Exception as exc:
            _calc_set_ops_status(
                f'Create folder failed: {_html.escape(str(exc))}',
                color='#d32f2f',
            )

    def calc_on_rename(button=None):
        try:
            src = _calc_collect_rename_source()
            src_resolved = _calc_resolve_within_root(src)
            current_before = _calc_current_dir().resolve()
            root = _calc_dir().resolve()
            new_name = _calc_clean_item_name(calc_rename_input.value)
            if src_resolved.name == new_name:
                raise ValueError('New name is identical to current name.')
            dst = _calc_resolve_within_root(src_resolved.parent / new_name)
            if dst.exists():
                raise ValueError(f'Target already exists: {new_name}')
            src_resolved.rename(dst)
            if state['current_path'] and current_before == src_resolved:
                state['current_path'] = dst.relative_to(root).as_posix()
            calc_rename_input.value = ''
            _calc_set_ops_status(
                f'Renamed to <code>{_html.escape(new_name)}</code>',
                color='#2e7d32',
            )
            calc_file_list.value = ()
            calc_file_list.value = ()
            _calc_refresh_related_explorers()
            return True
        except Exception as exc:
            _calc_set_ops_status(
                f'Rename failed: {_html.escape(str(exc))}',
                color='#d32f2f',
            )
            return False

    def calc_on_move_items(button=None):
        sources = _calc_collect_move_sources()
        if not sources:
            _calc_set_ops_status(
                'No source selected. Use Select mode or navigate into a folder.',
                color='#d32f2f',
            )
            return
        try:
            target_dir = _calc_parse_move_target_dir(calc_move_target_input.value)
        except Exception as exc:
            _calc_set_ops_status(
                f'Move failed: {_html.escape(str(exc))}',
                color='#d32f2f',
            )
            return

        current_before = _calc_current_dir().resolve()
        root = _calc_dir().resolve()
        moved, errors = [], []
        moved_current_to = None
        for src in sources:
            try:
                src_resolved = _calc_resolve_within_root(src)
                if src_resolved == target_dir or _calc_is_relative_to(target_dir, src_resolved):
                    raise ValueError('Cannot move a folder into itself.')
                dst = _calc_resolve_within_root(target_dir / src_resolved.name)
                if dst.exists():
                    raise ValueError(f'Target already exists: {src_resolved.name}')
                shutil.move(str(src_resolved), str(dst))
                moved.append(src_resolved.name)
                if state['current_path'] and current_before == src_resolved:
                    moved_current_to = dst
            except Exception as exc:
                errors.append(
                    f'{_html.escape(getattr(src, "name", str(src)))}: {_html.escape(str(exc))}'
                )

        if moved_current_to is not None:
            try:
                state['current_path'] = moved_current_to.resolve().relative_to(root).as_posix()
            except Exception:
                state['current_path'] = ''

        try:
            rel_target = target_dir.resolve().relative_to(root).as_posix()
            target_label = '/' if not rel_target else f'/{rel_target}'
        except Exception:
            target_label = str(target_dir)

        if moved and not errors:
            _calc_set_ops_status(
                f'Moved {len(moved)} item(s) to <code>{_html.escape(target_label)}</code>.',
                color='#2e7d32',
            )
        elif moved and errors:
            _calc_set_ops_status(
                (
                    f'Moved {len(moved)} item(s), failed: '
                    f'{"; ".join(errors)}'
                ),
                color='#e65100',
            )
        else:
            _calc_set_ops_status('; '.join(errors), color='#d32f2f')
        calc_move_target_input.value = ''
        calc_file_list.value = ()
        calc_file_list.value = ()
        _calc_refresh_related_explorers()

    def calc_on_explorer_new_folder(button=None):
        _calc_hide_rename_prompt()
        calc_new_folder_prompt_input.value = ''
        calc_new_folder_prompt_row.layout.display = 'flex'
        _calc_set_ops_status('Enter folder name and confirm with Create.', color='#555')

    def calc_on_new_folder_prompt_ok(button=None):
        raw_name = str(calc_new_folder_prompt_input.value or '').strip()
        calc_new_folder_prompt_row.layout.display = 'none'
        if not raw_name:
            _calc_set_ops_status('Folder name cannot be empty.', color='#d32f2f')
            return
        calc_new_folder_input.value = raw_name
        calc_on_new_folder()

    def calc_on_new_folder_prompt_cancel(button=None):
        _calc_hide_rename_prompt()
        _calc_set_ops_status('New folder canceled.', color='#555')

    def calc_on_explorer_start_rename(button=None):
        try:
            _calc_hide_rename_prompt()
            src = _calc_collect_rename_source()
            src_resolved = _calc_resolve_within_root(src)
            state['rename_source_path'] = str(src_resolved)
            calc_rename_prompt_input.value = src_resolved.name
            calc_rename_prompt_row.layout.display = 'flex'
            _calc_set_ops_status('Enter new name and confirm with OK.', color='#555')
        except Exception as exc:
            _calc_hide_rename_prompt()
            _calc_set_ops_status(
                f'Rename failed: {_html.escape(str(exc))}',
                color='#d32f2f',
            )

    def calc_on_explorer_start_duplicate(button=None):
        try:
            _calc_hide_rename_prompt()
            src_resolved = _calc_collect_duplicate_source()
            state['duplicate_source_path'] = str(src_resolved)
            calc_duplicate_prompt_input.value = src_resolved.name
            calc_duplicate_prompt_row.layout.display = 'flex'
            _calc_set_ops_status('Adjust folder name and confirm with Yes.', color='#555')
        except Exception as exc:
            _calc_hide_rename_prompt()
            _calc_set_ops_status(
                f'Duplicate failed: {_html.escape(str(exc))}',
                color='#d32f2f',
            )

    def calc_on_explorer_rename_ok(button=None):
        src_raw = str(state.get('rename_source_path') or '').strip()
        if not src_raw:
            _calc_set_ops_status('No item prepared for rename.', color='#d32f2f')
            return
        try:
            src_path = _calc_resolve_within_root(Path(src_raw))
            current_dir = _calc_current_dir().resolve()
            source_parent = src_path.parent.resolve()
            source_is_current_dir = src_path.resolve() == current_dir
            if source_parent != current_dir and not source_is_current_dir:
                raise ValueError('Open the source folder before confirming rename.')
            calc_action_source_input.value = src_path.name
            calc_rename_input.value = calc_rename_prompt_input.value
            success = calc_on_rename()
            if success:
                _calc_hide_rename_prompt()
        except Exception as exc:
            _calc_set_ops_status(
                f'Rename failed: {_html.escape(str(exc))}',
                color='#d32f2f',
            )

    def calc_on_explorer_duplicate_yes(button=None):
        src_raw = str(state.get('duplicate_source_path') or '').strip()
        if not src_raw:
            _calc_set_ops_status('No folder prepared for duplicate.', color='#d32f2f')
            return
        try:
            src_path = _calc_resolve_within_root(Path(src_raw))
            if not src_path.is_dir():
                raise ValueError('Duplicate works only for folders.')
            current_dir = _calc_current_dir().resolve()
            source_parent = src_path.parent.resolve()
            source_is_current_dir = src_path.resolve() == current_dir
            if source_parent != current_dir and not source_is_current_dir:
                raise ValueError('Open the source folder before confirming duplicate.')
            new_name = _calc_clean_item_name(calc_duplicate_prompt_input.value)
            dst = _calc_resolve_within_root(src_path.parent / new_name)
            if dst.exists():
                raise ValueError(f'Target already exists: {new_name}')
            shutil.copytree(str(src_path), str(dst))
            calc_file_list.value = ()
            calc_file_list.value = ()
            _calc_hide_rename_prompt()
            if source_is_current_dir:
                message = (
                    f'Duplicated to <code>{_html.escape(new_name)}</code>. '
                    'Go up one level to see it.'
                )
            else:
                message = f'Duplicated to <code>{_html.escape(new_name)}</code>.'
            _calc_set_ops_status(message, color='#2e7d32')
            _calc_refresh_related_explorers()
        except Exception as exc:
            _calc_set_ops_status(
                f'Duplicate failed: {_html.escape(str(exc))}',
                color='#d32f2f',
            )

    def calc_on_explorer_duplicate_cancel(button=None):
        _calc_hide_rename_prompt()
        _calc_set_ops_status('Duplicate canceled.', color='#555')

    def calc_on_explorer_rename_cancel(button=None):
        _calc_hide_rename_prompt()
        _calc_set_ops_status('Rename canceled.', color='#555')

    def _calc_clipboard_set(sources, mode):
        """Put *sources* on the shared clipboard with the given *mode* ('cut'/'copy')."""
        resolved_sources = []
        errors = []
        for src in sources:
            try:
                resolved_sources.append(src.resolve())
            except Exception as exc:
                errors.append(_html.escape(str(exc)))
        if not resolved_sources:
            _calc_set_ops_status(
                f'{mode.title()} failed: {"; ".join(errors) if errors else "No valid source selected."}',
                color='#d32f2f',
            )
            return
        normalized = [str(p) for p in resolved_sources]
        store = _calc_clipboard_store()
        store['paths'] = normalized
        store['mode'] = mode
        ctx.clipboard_paths = normalized
        ctx.clipboard_mode = mode
        _calc_update_explorer_action_state()
        verb = 'Cut' if mode == 'cut' else 'Copied'
        _calc_set_ops_status(
            f'{verb} {len(resolved_sources)} item(s). Open target folder and press Ctrl+V.',
            color='#1976d2',
        )

    def calc_on_explorer_cut(button=None):
        _calc_hide_rename_prompt()
        sources = _calc_collect_selected_sources_only()
        if not sources:
            _calc_set_ops_status('Select item(s) first.', color='#d32f2f')
            return
        _calc_clipboard_set(sources, 'cut')

    def calc_on_explorer_copy(button=None):
        _calc_hide_rename_prompt()
        sources = _calc_collect_selected_sources_only()
        if not sources:
            _calc_set_ops_status('Select item(s) first.', color='#d32f2f')
            return
        _calc_clipboard_set(sources, 'copy')

    def calc_on_explorer_paste(button=None):
        _calc_hide_rename_prompt()
        sources = _calc_clipboard_paths()
        if not sources:
            _calc_update_explorer_action_state()
            _calc_set_ops_status('Clipboard is empty. Use Ctrl+X or Ctrl+C first.', color='#d32f2f')
            return
        mode = _calc_clipboard_mode() or 'cut'
        target_dir = _calc_current_dir().resolve()
        target_dir.mkdir(parents=True, exist_ok=True)
        root = _calc_dir().resolve()
        done, errors = [], []
        remaining = []
        moved_current_to = None
        current_before = _calc_current_dir().resolve()
        for src in sources:
            try:
                src_resolved = src.resolve()
                if src_resolved == target_dir or _calc_is_relative_to(target_dir, src_resolved):
                    raise ValueError('Cannot paste a folder into itself.')
                dst = target_dir / src_resolved.name
                if dst.exists():
                    # Auto-number: name_2, name_3, ...
                    stem = dst.stem
                    suffix = dst.suffix
                    counter = 2
                    while dst.exists():
                        dst = target_dir / f'{stem}_{counter}{suffix}'
                        counter += 1
                if mode == 'copy':
                    if src_resolved.is_dir():
                        shutil.copytree(str(src_resolved), str(dst))
                    else:
                        shutil.copy2(str(src_resolved), str(dst))
                else:
                    shutil.move(str(src_resolved), str(dst))
                    if state['current_path'] and current_before == src_resolved:
                        moved_current_to = dst
                done.append(src_resolved.name)
            except Exception as exc:
                remaining.append(src)
                errors.append(
                    f'{_html.escape(getattr(src, "name", str(src)))}: {_html.escape(str(exc))}'
                )
        if mode == 'cut':
            normalized = [str(p) for p in remaining]
            store = _calc_clipboard_store()
            store['paths'] = normalized
            ctx.clipboard_paths = normalized
        # copy keeps clipboard intact
        _calc_update_explorer_action_state()

        if moved_current_to is not None:
            try:
                state['current_path'] = moved_current_to.resolve().relative_to(root).as_posix()
            except Exception:
                state['current_path'] = ''

        try:
            rel_target = target_dir.resolve().relative_to(root).as_posix()
            target_label = '/' if not rel_target else f'/{rel_target}'
        except Exception:
            target_label = str(target_dir)

        verb = 'Moved' if mode == 'cut' else 'Copied'
        if done and not errors:
            _calc_set_ops_status(
                f'{verb} {len(done)} item(s) to <code>{_html.escape(target_label)}</code>.',
                color='#2e7d32',
            )
        elif done and errors:
            _calc_set_ops_status(
                f'{verb} {len(done)} item(s), failed: {"; ".join(errors)}',
                color='#e65100',
            )
        else:
            _calc_set_ops_status('; '.join(errors), color='#d32f2f')
        calc_file_list.value = ()
        calc_file_list.value = ()
        _calc_refresh_related_explorers()

    # -- table extraction ---------------------------------------------------
    _CALC_TABLE_PRESET_SUFFIX = '.delfin-table.json'

    def _default_table_col_defs():
        return [{'name': 'Value 1', 'type': 'text', 'pattern': '', 'occ': 'last'}]

    def _normalize_table_preset_col_defs(columns):
        normalized = []
        for col in columns or []:
            if not isinstance(col, dict):
                continue
            col_type = str(col.get('type', 'text') or 'text').strip().lower()
            if col_type not in ('text', 'regex', 'json'):
                col_type = 'text'
            occ = str(col.get('occ', 'last') or 'last').strip().lower()
            if occ not in ('last', 'first', 'all'):
                occ = 'last'
            normalized.append({
                'name': str(col.get('name', '') or '').strip(),
                'type': col_type,
                'pattern': str(col.get('pattern', '') or '').strip(),
                'occ': occ,
            })
        return normalized or _default_table_col_defs()

    def _set_calc_table_preset_status(message='', color='#555'):
        if not message:
            calc_table_preset_status.value = ''
            return
        calc_table_preset_status.value = (
            f'<span style="color:{color};">{_html.escape(str(message))}</span>'
        )

    def _calc_is_table_preset_path(path):
        try:
            name = Path(path).name
        except Exception:
            name = str(path or '')
        return str(name or '').lower().endswith(_CALC_TABLE_PRESET_SUFFIX)

    def _normalize_calc_table_preset_name(raw_name=''):
        name = str(raw_name or '').strip()
        if not name:
            name = f'extract_table{_CALC_TABLE_PRESET_SUFFIX}'
        if '/' in name or '\\' in name:
            raise ValueError('Use only a file name, not a path.')
        if not name.lower().endswith(_CALC_TABLE_PRESET_SUFFIX):
            name += _CALC_TABLE_PRESET_SUFFIX
        return name

    def _table_preset_payload_from_ui():
        _collect_table_col_values()
        scope = str(calc_table_scope_dd.value or 'all').strip().lower()
        if scope not in ('all', 'selected'):
            scope = 'all'
        return {
            'kind': 'delfin_extract_table_preset',
            'version': 1,
            'file': str(calc_table_file_input.value or '').strip(),
            'scope': scope,
            'recursive': bool(calc_table_recursive_cb.value),
            'decimal_comma': bool(calc_table_decimal_comma_btn.value),
            'columns': _normalize_table_preset_col_defs(state.get('table_col_defs')),
        }

    def _apply_table_preset_payload(preset_payload, preset_name=''):
        payload = preset_payload if isinstance(preset_payload, dict) else {}
        scope = str(payload.get('scope', 'all') or 'all').strip().lower()
        if scope not in ('all', 'selected'):
            scope = 'all'
        calc_table_file_input.value = str(payload.get('file', '') or '').strip()
        calc_table_scope_dd.value = scope
        calc_table_recursive_cb.value = bool(payload.get('recursive', False))
        calc_table_decimal_comma_btn.value = bool(payload.get('decimal_comma', False))
        state['table_col_defs'] = _normalize_table_preset_col_defs(payload.get('columns'))
        _rebuild_table_col_rows()

        if preset_name:
            calc_table_preset_name.value = str(preset_name or '').strip()

    def _parse_calc_table_preset_text(raw_text, source_label=''):
        try:
            payload = json.loads(raw_text)
        except Exception as exc:
            raise ValueError(f'Invalid preset JSON: {exc}') from exc
        if not isinstance(payload, dict):
            raise ValueError('Preset file must contain a JSON object.')
        kind = str(payload.get('kind', '') or '').strip()
        if kind and kind != 'delfin_extract_table_preset':
            raise ValueError('Not a DELFIN extract-table preset file.')
        payload = dict(payload)
        payload['columns'] = _normalize_table_preset_col_defs(payload.get('columns'))
        payload['scope'] = str(payload.get('scope', 'all') or 'all').strip().lower()
        if payload['scope'] not in ('all', 'selected'):
            payload['scope'] = 'all'
        payload['recursive'] = bool(payload.get('recursive', False))
        payload['decimal_comma'] = bool(payload.get('decimal_comma', False))
        payload['file'] = str(payload.get('file', '') or '').strip()
        if not payload['file']:
            raise ValueError('Preset file is missing the table target file.')
        return payload

    def _open_calc_table_preset_file(full_path):
        raw_text = Path(full_path).read_text(encoding='utf-8')
        payload = _parse_calc_table_preset_text(raw_text, source_label=Path(full_path).name)
        calc_table_panel.layout.display = ''
        state['table_panel_active'] = True
        _apply_table_preset_payload(payload, preset_name=Path(full_path).name)
        _set_calc_table_preset_status(
            f'Loaded preset file {Path(full_path).name}. Click Run to extract the table.',
            color='#2e7d32',
        )
        calc_update_view()

    def _collect_table_col_values():
        for i, rw in enumerate(state['table_col_widgets']):
            if i < len(state['table_col_defs']):
                state['table_col_defs'][i]['name'] = rw['name'].value
                state['table_col_defs'][i]['type'] = rw['type_dd'].value
                state['table_col_defs'][i]['occ'] = rw['occ_dd'].value
                state['table_col_defs'][i]['pattern'] = rw['pattern'].value

    def _rebuild_table_col_rows():
        rows = []
        state['table_col_widgets'] = []
        for i, col in enumerate(state['table_col_defs']):
            name_w = widgets.Text(
                value=col.get('name', ''),
                placeholder='Column name',
                layout=widgets.Layout(width='100px', height='26px'),
            )
            type_dd = widgets.Dropdown(
                options=[('Text', 'text'), ('Regex', 'regex'), ('JSON path', 'json')],
                value=col.get('type', 'text'),
                layout=widgets.Layout(width='92px', height='26px'),
            )
            occ_dd = widgets.Dropdown(
                options=[('last', 'last'), ('first', 'first'), ('all', 'all')],
                value=col.get('occ', 'last'),
                layout=widgets.Layout(width='66px', height='26px'),
            )
            _pat_placeholders = {'regex': 'regex pattern', 'json': 'key.path', 'text': 'literal text'}
            pat_w = widgets.Text(
                value=col.get('pattern', ''),
                placeholder=_pat_placeholders.get(col.get('type', 'text'), 'literal text'),
                layout=widgets.Layout(flex='1 1 0', min_width='0', width='1px', height='26px'),
            )
            pat_w.add_class('calc-table-pattern-input')
            preset_dd = widgets.Dropdown(
                options=_tp_labels,
                value=_tp_labels[0],
                layout=widgets.Layout(width='180px', height='26px'),
            )
            rm_btn = widgets.Button(
                description='✕',
                layout=widgets.Layout(width='32px', height='26px'),
            )

            def _on_type_change(change, pw=pat_w):
                _ph = {'regex': 'regex pattern', 'json': 'key.path', 'text': 'literal text'}
                pw.placeholder = _ph.get(change['new'], 'literal text')

            def _on_preset(change, idx=i):
                label = change['new']
                if label == _tp_labels[0]:
                    return
                pi = _tp_labels.index(label)
                _collect_table_col_values()
                state['table_col_defs'][idx].update({
                    'name': _tp_names[pi],
                    'type': _tp_types[pi],
                    'pattern': _tp_pats[pi],
                    'occ': _tp_occs[pi],
                })
                if _tp_types[pi] == 'json' and not calc_table_file_input.value.strip():
                    calc_table_file_input.value = 'DELFIN_Data.json'
                _rebuild_table_col_rows()

            def _on_remove(b, idx=i):
                _collect_table_col_values()
                if len(state['table_col_defs']) > 1:
                    state['table_col_defs'].pop(idx)
                _rebuild_table_col_rows()

            type_dd.observe(_on_type_change, names='value')
            preset_dd.observe(_on_preset, names='value')
            rm_btn.on_click(_on_remove)

            row_box = widgets.HBox(
                [name_w, type_dd, occ_dd, pat_w, preset_dd, rm_btn],
                layout=widgets.Layout(gap='4px', align_items='center', width='100%'),
            )
            row_box.add_class('calc-table-col-row')
            rows.append(row_box)
            state['table_col_widgets'].append({
                'name': name_w, 'type_dd': type_dd, 'occ_dd': occ_dd,
                'pattern': pat_w, 'preset_dd': preset_dd,
            })
        calc_table_cols_box.children = tuple(rows)

    def _json_collect_key_values(node, key):
        matches = []
        stack = [node]
        while stack:
            current = stack.pop()
            if isinstance(current, dict):
                if key in current:
                    matches.append(current[key])
                values = list(current.values())
                for value in reversed(values):
                    stack.append(value)
            elif isinstance(current, list):
                for item in reversed(current):
                    stack.append(item)
        return matches

    def _normalize_json_path(path):
        """Normalize common JSONPath-like prefixes to the local dot syntax."""
        normalized = str(path or '').strip()
        if not normalized:
            return ''
        if normalized.startswith('$'):
            normalized = normalized[1:]
        return normalized.lstrip('.')

    def _json_split_path(path):
        """Split JSON path by dots while preserving bracket expressions."""
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
                if ch == '\\':
                    escaped = True
                elif ch == quote:
                    quote = None
                continue
            if ch in ("'", '"') and bracket_depth > 0:
                quote = ch
                buf.append(ch)
                continue
            if ch == '[':
                bracket_depth += 1
                buf.append(ch)
                continue
            if ch == ']':
                if bracket_depth > 0:
                    bracket_depth -= 1
                buf.append(ch)
                continue
            if ch == '.' and bracket_depth == 0:
                token = ''.join(buf).strip()
                if token:
                    parts.append(token)
                buf = []
                continue
            buf.append(ch)
        tail = ''.join(buf).strip()
        if tail:
            parts.append(tail)
        return parts

    def _json_parse_value(raw):
        token = str(raw).strip()
        if len(token) >= 2 and token[0] == token[-1] and token[0] in ('"', "'"):
            return token[1:-1]
        low = token.lower()
        if low == 'true':
            return True
        if low == 'false':
            return False
        if low == 'null':
            return None
        try:
            if any(ch in token for ch in ('.', 'e', 'E')):
                return float(token)
            return int(token)
        except Exception:
            return token

    def _json_compare_values(lhs, op, rhs):
        if op in ('==', '!='):
            equal = (lhs == rhs)
            if not equal:
                try:
                    equal = (float(lhs) == float(rhs))
                except Exception:
                    equal = (str(lhs) == str(rhs))
            return equal if op == '==' else (not equal)
        try:
            lval = float(lhs)
            rval = float(rhs)
        except Exception:
            lval = str(lhs)
            rval = str(rhs)
        if op == '>':
            return lval > rval
        if op == '>=':
            return lval >= rval
        if op == '<':
            return lval < rval
        if op == '<=':
            return lval <= rval
        return False

    def _json_item_matches_filter(item, expression):
        if not isinstance(item, dict):
            return False
        clauses = [
            c.strip()
            for c in re.split(r'\s*(?:&&|\band\b)\s*', expression)
            if c.strip()
        ]
        if not clauses:
            return False
        for clause in clauses:
            m = re.match(r'^@\.(\w+)\s*(==|!=|>=|<=|>|<)\s*(.+)$', clause)
            if not m:
                return False
            key, op, raw_rhs = m.groups()
            rhs = _json_parse_value(raw_rhs)
            lhs = item.get(key)
            if not _json_compare_values(lhs, op, rhs):
                return False
        return True

    def _json_parse_part(part):
        """
        Parse one path part and return (key, selector_kind, selector_value).
        Supported selectors: [*], [index], [?(...)].
        """
        text = part.strip()
        if not text:
            return '', None, None
        if '[' not in text or not text.endswith(']'):
            return text, None, None
        bracket_pos = text.find('[')
        key = text[:bracket_pos].strip()
        selector = text[bracket_pos + 1:-1].strip()
        if selector == '*':
            return key, 'wildcard', None
        if re.fullmatch(r'-?\d+', selector):
            return key, 'index', int(selector)
        if selector.startswith('?(') and selector.endswith(')'):
            return key, 'filter', selector[2:-1].strip()
        return key, None, None

    def _json_apply_selector(values, selector_kind, selector_value):
        if selector_kind is None:
            return values
        selected = []
        for value in values:
            if selector_kind == 'wildcard':
                if isinstance(value, list):
                    selected.extend(value)
                elif isinstance(value, dict):
                    selected.extend(value.values())
            elif selector_kind == 'index':
                if isinstance(value, list):
                    try:
                        selected.append(value[selector_value])
                    except Exception:
                        pass
            elif selector_kind == 'filter':
                if isinstance(value, list):
                    for item in value:
                        if _json_item_matches_filter(item, selector_value):
                            selected.append(item)
                elif _json_item_matches_filter(value, selector_value):
                    selected.append(value)
        return selected

    def _json_step(node, part):
        key, selector_kind, selector_value = _json_parse_part(part)
        lower_key = key.lower()
        if isinstance(node, dict):
            if key == '*':
                matches = list(node.values())
            elif key in node:
                matches = [node[key]]
            else:
                # Convenience fallback: search this key recursively.
                matches = _json_collect_key_values(node, key)
            return _json_apply_selector(matches, selector_kind, selector_value)
        if isinstance(node, list):
            if key in ('*', '[]'):
                matches = list(node)
            elif lower_key == 'last':
                matches = [node[-1]] if node else []
            elif lower_key == 'first':
                matches = [node[0]] if node else []
            elif key == '':
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
        typ = col_def.get('type', 'regex')
        pat = col_def.get('pattern', '').strip()
        occ = col_def.get('occ', 'last')
        if not pat:
            return ['—']
        try:
            if typ in ('regex', 'text'):
                rx = re.compile(re.escape(pat) if typ == 'text' else pat)
                if rx.groups > 0:
                    # explicit capture group — no line tracking possible
                    raw_matches = rx.finditer(text)
                    values = []  # list of (value_str, line_num_or_None)
                    for m in raw_matches:
                        g = m.group(1)
                        val = str(g).strip() if g is not None else '—'
                        line_num = text[:m.start()].count('\n') + 1
                        values.append((val, line_num))
                    if not values:
                        return ['—']
                else:
                    values = []  # list of (value_str, line_num)
                    for m in rx.finditer(text):
                        line_num = text[:m.start()].count('\n') + 1
                        line_end = text.find('\n', m.end())
                        if line_end == -1:
                            line_end = len(text)
                        tail = text[m.end():line_end]
                        nums = [nm.group(0) for nm in _table_number_re.finditer(tail)]
                        if nums:
                            values.append((_normalize_table_number(nums[0]), line_num))
                            continue
                        # wrapped lines: inspect a short window after the match
                        near = text[m.end():min(len(text), m.end() + 240)]
                        near_nums = [nm.group(0) for nm in _table_number_re.finditer(near)]
                        if near_nums:
                            values.append((_normalize_table_number(near_nums[0]), line_num))
                            continue
                        # fallback: keep old behavior if no number follows
                        values.append((m.group(0).strip(), line_num))
                    if not values:
                        return ['—']
                if occ == 'all':
                    return values  # list of (value_str, line_num)
                if occ == 'first':
                    return [values[0][0]]
                return [values[-1][0]]
            elif typ == 'json':
                if json_data is None:
                    return ['—']
                values = _json_extract_path_values(json_data, pat)
                if not values:
                    return ['—']
                rendered = [str(v) for v in values]
                if occ == 'all':
                    return rendered
                if occ == 'first':
                    return [rendered[0]]
                return [rendered[-1]]
        except Exception as exc:
            return [f'err:{exc}']
        return ['—']

    def _format_table_output_value(raw_value):
        value = str(raw_value)
        if not calc_table_decimal_comma_btn.value:
            return value
        if not _table_number_full_re.fullmatch(value.strip()):
            return value
        return value.replace('.', ',')

    def _render_extract_table_html(headers, rows):
        th = ''.join(
            f'<th style="padding:4px 10px;border:1px solid #ddd;'
            f'background:#f0f4f8;white-space:nowrap;text-align:left">'
            f'{_html.escape(str(h))}</th>'
            for h in headers
        )
        trs = ''
        for i, row in enumerate(rows):
            bg = '#ffffff' if i % 2 == 0 else '#f7f9fc'
            tds = ''.join(
                f'<td style="padding:4px 10px;border:1px solid #ddd;'
                f'white-space:nowrap">{_html.escape(str(v))}</td>'
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

    def calc_on_table_toggle(button):
        if calc_table_panel.layout.display == 'none':
            calc_table_panel.layout.display = ''
            state['table_panel_active'] = True
            if not str(calc_table_preset_name.value or '').strip():
                calc_table_preset_name.value = f'extract_table{_CALC_TABLE_PRESET_SUFFIX}'
            _rebuild_table_col_rows()
        else:
            calc_table_panel.layout.display = 'none'
            state['table_panel_active'] = False
        calc_update_view()

    def calc_on_table_close(button):
        calc_table_panel.layout.display = 'none'
        state['table_panel_active'] = False
        calc_update_view()

    def calc_on_table_add_col(button):
        _collect_table_col_values()
        n = len(state['table_col_defs']) + 1
        state['table_col_defs'].append(
            {'name': f'Value {n}', 'type': 'regex', 'pattern': '', 'occ': 'last'}
        )
        _rebuild_table_col_rows()

    def calc_on_table_preset_save(button=None):
        try:
            preset_name = _normalize_calc_table_preset_name(calc_table_preset_name.value)
        except Exception as exc:
            _set_calc_table_preset_status(str(exc), color='#d32f2f')
            return
        base_dir = (
            _calc_dir() / state['current_path']
            if state['current_path'] else _calc_dir()
        )
        preset_path = base_dir / preset_name
        try:
            preset_path.write_text(
                json.dumps(_table_preset_payload_from_ui(), indent=2, ensure_ascii=True) + '\n',
                encoding='utf-8',
            )
        except Exception as exc:
            _set_calc_table_preset_status(f'Could not save preset file: {exc}', color='#d32f2f')
            return
        calc_table_preset_name.value = preset_name
        calc_list_directory()
        _set_calc_table_preset_status(
            f'Saved preset file {preset_name} in {base_dir}.',
            color='#2e7d32',
        )

    def calc_on_table_run(button):
        _collect_table_col_values()
        target_file = calc_table_file_input.value.strip()
        if not target_file:
            calc_table_status.value = '<span style="color:#d32f2f;">Please enter a filename.</span>'
            return
        base_dir = (
            _calc_dir() / state['current_path']
            if state['current_path'] else _calc_dir()
        )
        scope = calc_table_scope_dd.value
        if scope == 'selected':
            real = _calc_selected_labels()
            if not real:
                calc_table_status.value = (
                    '<span style="color:#d32f2f;">No folders selected.</span>'
                )
                return
            folders = []
            for label in real:
                p = base_dir / _calc_label_to_name(label)
                if p.is_dir():
                    folders.append(p)
        else:
            try:
                folders = sorted(
                    [e for e in base_dir.iterdir() if e.is_dir()],
                    key=lambda x: x.name.lower(),
                )
            except Exception:
                folders = []
        if not folders:
            calc_table_status.value = '<span style="color:#d32f2f;">No folders found.</span>'
            return
        cols = state['table_col_defs']
        headers = ['Job']
        for i, c in enumerate(cols):
            headers.append(c.get('name', f'Col {i+1}'))
            if c.get('occ') == 'all':
                headers.append(c.get('name', f'Col {i+1}') + ' (Line)')
        _n_header_cols = len(headers) - 1  # excluding Job
        rows = []
        missing = 0
        for folder in folders:
            direct = folder / target_file
            if direct.is_file():
                file_paths = [direct]
            elif calc_table_recursive_cb.value:
                try:
                    file_paths = sorted(folder.rglob(target_file))
                except Exception:
                    file_paths = []
            else:
                file_paths = []
            if not file_paths:
                rows.append([folder.name] + ['—'] * _n_header_cols)
                missing += 1
                continue
            for file_path in file_paths:
                try:
                    job_label = str(file_path.parent.relative_to(base_dir))
                except ValueError:
                    job_label = folder.name
                try:
                    content = file_path.read_text(errors='replace')
                except Exception:
                    rows.append([job_label] + ['err'] * _n_header_cols)
                    continue
                json_data = None
                if file_path.suffix.lower() == '.json':
                    try:
                        json_data = json.loads(content)
                    except Exception:
                        pass
                col_values = [_extract_values(c, content, json_data) for c in cols]
                row_count = max((len(vs) for vs in col_values), default=1)
                for row_idx in range(row_count):
                    row = [job_label]
                    for c, values in zip(cols, col_values):
                        is_all = c.get('occ') == 'all'
                        if not values or values == ['—']:
                            row.append('—')
                            if is_all:
                                row.append('—')
                        elif row_idx < len(values):
                            v = values[row_idx]
                            if isinstance(v, tuple):
                                row.append(_format_table_output_value(v[0]))
                                if is_all:
                                    row.append(str(v[1]) if v[1] is not None else '—')
                            else:
                                row.append(_format_table_output_value(v))
                                if is_all:
                                    row.append('—')
                        else:
                            row.append('—')
                            if is_all:
                                row.append('—')
                    rows.append(row)
        calc_table_output.value = _render_extract_table_html(headers, rows)
        import io, csv as _csv
        buf = io.StringIO()
        csv_rows = [['' if c == '—' else c for c in row] for row in rows]
        _csv.writer(buf, delimiter=';').writerows([headers] + csv_rows)
        state['table_csv_data'] = buf.getvalue()
        calc_table_csv_btn.layout.display = 'inline-flex'
        folder_count = len(folders)
        row_count = len(rows)
        msg = f'<span style="color:#2e7d32;">{folder_count} folder(s) processed'
        if row_count != folder_count:
            msg += f', {row_count} row(s) generated'
        if missing:
            msg += f', {missing} without {_html.escape(target_file)}'
        calc_table_status.value = msg + '.</span>'

    def calc_on_table_csv_download(button):
        data = state.get('table_csv_data', '')
        if not data:
            return
        import base64 as _b64
        b64 = _b64.b64encode(data.encode('utf-8-sig')).decode()
        _js = (
            'var a=document.createElement("a");'
            f'a.href="data:text/csv;charset=utf-8;base64,{b64}";'
            'a.download="extract_table.csv";'
            'document.body.appendChild(a);a.click();document.body.removeChild(a);'
        )
        from IPython.display import Javascript as _JS
        with ctx.js_output:
            display(_JS(_js))

    def calc_on_sort_change(change):
        saved_filter = calc_folder_search.value
        calc_list_directory()
        if saved_filter:
            calc_folder_search.value = saved_filter
            calc_filter_file_list()

    # -- item open logic (shared by dblclick and single-click on files) ------
    def _calc_open_item(selected):
        """Open/display a single item given its label string."""
        if not selected or selected.startswith('('):
            return

        name = _calc_label_to_name(selected)
        full_path = (
            (_calc_dir() / state['current_path'] / name)
            if state['current_path']
            else (_calc_dir() / name)
        )

        # Directories are only navigated via dblclick, not here.
        if full_path.is_dir():
            return

        if _is_archive_tab and _calc_is_table_preset_path(full_path):
            try:
                _open_calc_table_preset_file(full_path)
            except Exception as exc:
                _set_calc_table_preset_status(f'Could not load preset file: {exc}', color='#d32f2f')
                calc_table_panel.layout.display = ''
                state['table_panel_active'] = True
                calc_update_view()
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
            _calc_clear_main_viewer_state(reset_view_state=False)
            calc_mol_container.layout.display = 'none'
        calc_content_area.layout.display = 'none'
        calc_content_label.layout.display = 'none'
        if not keep_previous_viewer_during_load:
            _set_view_toggle(False, True)
        calc_file_info.value = ''
        calc_path_display.value = ''
        calc_download_status.value = ''
        calc_xyz_batch_status.value = ''
        state['xyz_batch_marked_paths'] = []
        calc_print_mode_status.value = ''
        calc_mo_plot_status.value = ''
        state['file_content'] = ''
        state['selected_file_path'] = None
        state['selected_file_size'] = 0
        state['file_is_preview'] = False
        state['file_preview_note'] = ''
        state['file_chunk_start'] = 0
        state['file_chunk_end'] = 0
        state['chunk_dom_initialized'] = False
        _calc_stop_xyz_playback(update_button=True)
        _calc_hide_chunk_controls()
        calc_reset_recalc_state()
        calc_reset_xyz_workflow_state()
        _calc_preselect_show(False)
        _calc_show_xyz_batch_panel(False)
        _calc_show_print_mode_panel(False)
        _calc_show_mo_plot_panel(False)

        # Reset xyz trajectory controls
        state['xyz_frames'].clear()
        state['xyz_current_frame'][0] = 0
        state['traj_viewer_ready'] = False
        calc_xyz_controls.layout.display = 'none'
        calc_coord_controls.layout.display = 'none'
        _calc_update_xyz_traj_control_state()
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
                    _calc_update_xyz_traj_control_state()
                    _render_3dmol(xyz_data)
                    _calc_apply_traj_style()
                else:
                    calc_file_info.value = f'<b>{_html.escape(name)}</b> (could not parse)'
            except Exception as e:
                calc_file_info.value = f'<b>{_html.escape(name)}</b> (error: {e})'
            return

        # --- XYZ file (with trajectory support) ---
        if suffix == '.xyz':
            if size > CALC_XYZ_MAX_READ_BYTES:
                # Large trajectory: read only enough for the first frame(s)
                # to show in 3D viewer, and use chunked text for the rest.
                _set_view_toggle(True, False)
                calc_mol_container.layout.display = 'block'
                calc_content_area.layout.display = 'none'
                calc_update_view()
                try:
                    head_bytes = min(CALC_TEXT_CHUNK_BYTES, size)
                    with open(str(full_path), 'r', errors='ignore') as fh:
                        head_text = fh.read(head_bytes)
                    head_frames = parse_xyz_frames(head_text)
                    _calc_load_text_preview_chunk(full_path, size, 0)
                    calc_copy_btn.disabled = False
                    calc_copy_path_btn.disabled = False
                    if head_frames:
                        state['xyz_frames'].clear()
                        state['xyz_frames'].extend(head_frames)
                        state['xyz_current_frame'][0] = 0
                        n_head = len(head_frames)
                        calc_xyz_frame_input.max = n_head
                        calc_xyz_frame_input.value = 1
                        calc_xyz_frame_total.value = f"<b>/ {n_head}+</b>"
                        calc_xyz_controls.layout.display = 'flex'
                        calc_coord_controls.layout.display = 'none'
                        _calc_update_xyz_traj_control_state()
                        calc_update_xyz_viewer(initial_load=True)
                        calc_file_info.value = (
                            f'<b><span style="word-break:break-all;">{_html.escape(name)}</span></b>'
                            f' ({size_str}, showing first {n_head} frames, file too large for full load)'
                        )
                    else:
                        calc_file_info.value = (
                            f'<b><span style="word-break:break-all;">{_html.escape(name)}</span></b>'
                            f' ({size_str}, chunked view)'
                        )
                    _set_view_toggle(False, False)
                except Exception as e:
                    calc_set_message(f'Error: {e}')
                return
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
                        _calc_update_xyz_traj_control_state()
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
                    _calc_update_xyz_traj_control_state()
                    _calc_set_rmsd_available(False)
                    calc_update_xyz_viewer(initial_load=True)
                else:
                    calc_file_info.value = (
                        f'<b><span style="word-break:break-all;">{_html.escape(name)}</span></b>'
                        f' ({size_str}, {len(lines)} lines)'
                    )
                    if state['xyz_frames']:
                        comment = state['xyz_frames'][0][0]
                        state['xyz_current_frame'][0] = 0
                        calc_xyz_frame_input.max = len(state['xyz_frames'])
                        calc_xyz_frame_input.value = 1
                        calc_xyz_frame_total.value = f"<b>/ {len(state['xyz_frames'])}</b>"
                        calc_xyz_controls.layout.display = 'flex'
                        calc_coord_controls.layout.display = 'none'
                        calc_xyz_frame_label.value = format_xyz_comment_label(
                            comment,
                            frames=state['xyz_frames'],
                        )
                        _calc_update_xyz_traj_control_state()
                        _calc_set_rmsd_available(True)
                    else:
                        calc_xyz_controls.layout.display = 'none'
                        _calc_update_xyz_traj_control_state()
                        _calc_set_rmsd_available(False)
                        calc_xyz_frame_label.value = ''
                    _render_3dmol(content)
                    _calc_apply_traj_style()
                calc_render_content(scroll_to='top')
            except Exception as e:
                calc_set_message(f'Error: {e}')
            return

        # --- PNG image ---
        if suffix == '.png':
            _calc_set_png_button_mode(main=False)
            if size > CALC_IMAGE_MAX_READ_BYTES:
                calc_file_info.value = (
                    f'<b><span style="word-break:break-all;">{_html.escape(name)}</span></b>'
                    f' ({size_str})'
                )
                calc_set_message(f'Image too large for inline display ({size_str}).')
                return
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
            if size > CALC_CUBE_MAX_READ_BYTES:
                calc_file_info.value = (
                    f'<b><span style="word-break:break-all;">{_html.escape(name)}</span></b>'
                    f' ({size_str})'
                )
                calc_set_message(
                    f'Cube file too large for live preview ({size_str}).\n\n'
                    f'Use chunked text view instead.'
                )
                _calc_load_text_preview_chunk(full_path, size, 0)
                calc_copy_btn.disabled = False
                calc_copy_path_btn.disabled = False
                calc_update_view()
                return
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
                    _calc_set_png_button_mode(main=False)
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
        _BINARY_SUFFIXES = {
            '.gbw', '.cis', '.densities', '.tmp', '.rwf', '.chk', '.fchk',
            '.wfn', '.wfx', '.molden', '.nat', '.hess', '.engrad',
            '.zip', '.gz', '.tar', '.bz2', '.xz', '.7z',
            '.so', '.o', '.a', '.exe', '.dll', '.pyc', '.pyo',
            '.db', '.sqlite', '.sqlite3',
        }
        if suffix in _BINARY_SUFFIXES:
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
                elif suffix in ['.out', '.log']:
                    try:
                        if size > CALC_LOG_XYZ_EXTRACT_MAX:
                            # Read only the tail of large logs (XYZ blocks are near the end)
                            with open(str(full_path), 'r', errors='ignore') as fh:
                                fh.seek(max(0, size - CALC_TEXT_CHUNK_BYTES))
                                xyz_source = fh.read()
                        else:
                            xyz_source = full_path.read_text(errors='ignore')
                        xyz_content = calc_extract_orca_xyz_block(xyz_source)
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
                                _calc_update_xyz_traj_control_state()
                                if len(state['xyz_frames']) > 1:
                                    calc_update_xyz_viewer(initial_load=True)
                                else:
                                    calc_xyz_frame_label.value = ''
                                    _render_3dmol(xyz_content)
                                    _calc_apply_traj_style()
                            else:
                                calc_xyz_controls.layout.display = 'none'
                                _calc_update_xyz_traj_control_state()
                            _set_view_toggle(False, False)
                    except Exception:
                        pass
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
                            _calc_update_xyz_traj_control_state()
                        else:
                            calc_xyz_controls.layout.display = 'none'
                            _calc_update_xyz_traj_control_state()
                        calc_xyz_frame_label.value = ''
                        _render_3dmol(xyz_content)
                        _calc_apply_traj_style()
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
                            _calc_update_xyz_traj_control_state()
                        else:
                            calc_xyz_controls.layout.display = 'none'
                            _calc_update_xyz_traj_control_state()
                        calc_xyz_frame_label.value = ''
                        _render_3dmol(xyz_content)
                        _calc_apply_traj_style()
                        _set_view_toggle(False, False)
                    else:
                        _set_view_toggle(False, True)
                except Exception as exc:
                    _set_view_toggle(False, True)
                    calc_set_message(f'Error: {exc}')
                calc_mol_container.layout.display = 'none'
                calc_content_area.layout.display = 'block'
            elif suffix in ['.out', '.log']:
                try:
                    xyz_content = calc_extract_orca_xyz_block(content)
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
                            _calc_update_xyz_traj_control_state()
                            if len(state['xyz_frames']) > 1:
                                calc_update_xyz_viewer(initial_load=True)
                            else:
                                calc_xyz_frame_label.value = ''
                                _render_3dmol(xyz_content)
                                _calc_apply_traj_style()
                        else:
                            calc_xyz_controls.layout.display = 'none'
                            _calc_update_xyz_traj_control_state()
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

    # -- double-click handler (navigate into dir, or open file) ---------------
    def calc_on_dblclick(change):
        raw = (change['new'] or '').strip()
        calc_dblclick_input.value = ''
        chosen_label, full_path = _calc_resolve_active_label(raw)
        if not chosen_label or full_path is None:
            return
        name = _calc_label_to_name(chosen_label)
        if full_path.is_dir():
            _calc_hide_rename_prompt()
            state['current_path'] = (
                f'{state["current_path"]}/{name}' if state['current_path'] else name
            )
            calc_list_directory()
            calc_file_list.value = ()
            return
        if full_path.is_file():
            _calc_open_item(chosen_label)

    # -- selection-change handler (show file content on single-click) --------
    def calc_on_selection_change(change):
        _calc_update_explorer_action_state()
        labels = [label for label in (change['new'] or ()) if label and not label.startswith('(')]
        if len(labels) == 1:
            name = _calc_label_to_name(labels[0])
            full_path = (
                (_calc_dir() / state['current_path'] / name)
                if state['current_path']
                else (_calc_dir() / name)
            )
            if full_path.is_file():
                _calc_open_item(labels[0])

    # -- keyboard action handler (from JS bridge) ----------------------------
    def calc_on_keyboard_action(change):
        action = (change['new'] or '').strip()
        calc_keyboard_action_input.value = ''
        if not action:
            return
        if action == 'copy':
            calc_on_explorer_copy()
        elif action == 'cut':
            calc_on_explorer_cut()
        elif action == 'paste':
            calc_on_explorer_paste()
        elif action == 'delete':
            calc_on_delete_click()
        elif action == 'rename':
            calc_on_explorer_start_rename()
        elif action == 'open':
            labels = _calc_selected_labels()
            if len(labels) == 1:
                calc_on_dblclick({'new': labels[0]})
        elif action == 'deselect':
            calc_file_list.value = ()

    # -- wiring -------------------------------------------------------------
    calc_back_btn.on_click(calc_on_back)
    calc_home_btn.on_click(calc_on_home)
    calc_refresh_btn.on_click(calc_on_refresh)
    calc_path_input.observe(calc_on_path_change, names='value')
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
    calc_submit_xyz_workflow_btn.on_click(_calc_submit_xyz_browser_workflow)
    calc_delete_btn.on_click(calc_on_delete_click)
    calc_delete_yes_btn.on_click(calc_on_delete_yes)
    calc_delete_no_btn.on_click(calc_on_delete_no)
    calc_move_archive_btn.on_click(calc_on_move_archive_click)
    calc_back_to_calculations_btn.on_click(calc_on_back_to_calculations)
    calc_ssh_transfer_btn.on_click(calc_on_ssh_transfer_click)
    calc_transfer_jobs_btn.on_click(calc_on_transfer_jobs_toggle)
    calc_move_archive_yes_btn.on_click(calc_on_move_archive_yes)
    calc_move_archive_no_btn.on_click(calc_on_move_archive_no)
    calc_transfer_jobs_refresh_btn.on_click(calc_on_transfer_jobs_refresh)
    calc_explorer_new_btn.on_click(calc_on_explorer_new_folder)
    calc_explorer_rename_btn.on_click(calc_on_explorer_start_rename)
    calc_duplicate_btn.on_click(calc_on_explorer_start_duplicate)
    calc_rename_prompt_ok_btn.on_click(calc_on_explorer_rename_ok)
    calc_rename_prompt_cancel_btn.on_click(calc_on_explorer_rename_cancel)
    calc_duplicate_prompt_yes_btn.on_click(calc_on_explorer_duplicate_yes)
    calc_duplicate_prompt_cancel_btn.on_click(calc_on_explorer_duplicate_cancel)
    calc_new_folder_prompt_ok_btn.on_click(calc_on_new_folder_prompt_ok)
    calc_new_folder_prompt_cancel_btn.on_click(calc_on_new_folder_prompt_cancel)
    calc_new_folder_btn.on_click(calc_on_new_folder)
    calc_rename_btn.on_click(calc_on_rename)
    calc_move_btn.on_click(calc_on_move_items)
    calc_hidden_upload.observe(calc_on_hidden_upload, names='value')
    calc_upload_seq_input.observe(calc_on_upload_bridge_seq, names='value')
    calc_upload_trigger_btn.on_click(calc_on_upload_trigger_click)
    for _txt, _handler in (
        (calc_new_folder_input, calc_on_new_folder),
        (calc_rename_input, calc_on_rename),
        (calc_move_target_input, calc_on_move_items),
    ):
        if hasattr(_txt, 'on_submit'):
            try:
                _txt.on_submit(_handler)
            except Exception:
                pass
    if hasattr(calc_rename_prompt_input, 'on_submit'):
        try:
            calc_rename_prompt_input.on_submit(calc_on_explorer_rename_ok)
        except Exception:
            pass
    if hasattr(calc_duplicate_prompt_input, 'on_submit'):
        try:
            calc_duplicate_prompt_input.on_submit(calc_on_explorer_duplicate_yes)
        except Exception:
            pass
    if hasattr(calc_new_folder_prompt_input, 'on_submit'):
        try:
            calc_new_folder_prompt_input.on_submit(calc_on_new_folder_prompt_ok)
        except Exception:
            pass
    if _is_archive_tab:
        calc_table_btn.on_click(calc_on_table_toggle)
        calc_table_close_btn.on_click(calc_on_table_close)
        calc_table_add_col_btn.on_click(calc_on_table_add_col)
        calc_table_run_btn.on_click(calc_on_table_run)
        calc_table_csv_btn.on_click(calc_on_table_csv_download)
        calc_table_preset_save_btn.on_click(calc_on_table_preset_save)
    calc_file_list.observe(calc_on_selection_change, names='value')
    calc_dblclick_input.observe(calc_on_dblclick, names='value')
    calc_xyz_batch_dblclick_input.observe(calc_on_xyz_batch_dblclick, names='value')
    calc_xyz_batch_toggle_input.observe(calc_on_xyz_batch_toggle, names='value')
    calc_keyboard_action_input.observe(calc_on_keyboard_action, names='value')
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
    calc_xyz_loop_checkbox.observe(calc_on_xyz_loop_change, names='value')
    calc_xyz_fps_input.observe(calc_on_xyz_fps_change, names='value')
    calc_xyz_play_btn.observe(calc_on_xyz_play_change, names='value')
    calc_xyz_copy_btn.on_click(calc_on_xyz_copy)
    calc_xyz_png_btn.on_click(calc_on_view_png)
    calc_rmsd_run_btn.on_click(calc_on_rmsd_run)
    calc_rmsd_hide_btn.on_click(calc_on_rmsd_hide)
    calc_coord_copy_btn.on_click(calc_on_coord_copy)
    calc_coord_png_btn.on_click(calc_on_view_png)
    calc_view_png_btn.on_click(calc_on_view_png)
    calc_copy_btn.on_click(calc_on_content_copy)
    calc_copy_path_btn.on_click(calc_on_path_copy)
    calc_download_btn.on_click(calc_on_download)
    calc_xyz_batch_mark_btn.on_click(calc_on_xyz_batch_mark)
    calc_xyz_batch_select_all_btn.on_click(calc_on_xyz_batch_select_all)
    calc_xyz_batch_up_btn.on_click(calc_on_xyz_batch_up)
    calc_xyz_batch_root_btn.on_click(calc_on_xyz_batch_root)
    calc_xyz_batch_refresh_btn.on_click(calc_on_xyz_batch_refresh)
    calc_xyz_batch_build_btn.on_click(calc_on_xyz_batch_build)
    calc_xyz_batch_copy_btn.on_click(calc_on_xyz_batch_copy)
    calc_print_mode_plot_btn.on_click(calc_on_print_mode_plot)
    calc_mo_plot_alpha_cb.observe(_calc_update_mo_spin_button_styles, names='value')
    calc_mo_plot_beta_cb.observe(_calc_update_mo_spin_button_styles, names='value')
    _calc_update_mo_spin_button_styles()
    calc_mo_plot_btn.on_click(calc_on_mo_plot)
    calc_nmr_submit_btn.on_click(_calc_on_nmr_submit)
    calc_censo_nmr_submit_btn.on_click(_calc_on_censo_nmr_submit)
    calc_report_btn.on_click(calc_on_report)
    disable_spellcheck(ctx, class_name='calc-search-input')
    disable_spellcheck(ctx, class_name='calc-table-preset-name')
    _calc_update_xyz_traj_control_state()

    # -- initialise ---------------------------------------------------------
    _calc_register_explorer_refresh()
    if _remote_archive_enabled:
        _calc_render_transfer_jobs()

    if ctx.calc_dir.exists():
        calc_list_directory()
    else:
        calc_file_list.options = ['(calc folder not found)']

    _explorer_interactions_js = r"""
    (function(){
        if (window._delfinExplorerInteractionsReady) return;
        window._delfinExplorerInteractionsReady = true;
        console.log('[DELFIN] Explorer interactions JS loaded');

        function _labelText(opt){
            if (!opt) return '';
            return String(opt.textContent || opt.innerText || '');
        }
        function _labelToName(text){
            var t = String(text || '').trim();
            var idx = t.indexOf(' ');
            return idx >= 0 ? t.slice(idx + 1).trim() : t;
        }
        function _isFolderLabel(text){
            var t = String(text || '');
            return t.indexOf('📂 ') === 0 || t.indexOf('✅ ') === 0;
        }
        function _eventPointNode(e){
            if (document.elementFromPoint && e && typeof e.clientX === 'number' && typeof e.clientY === 'number') {
                var pointed = document.elementFromPoint(e.clientX, e.clientY);
                if (pointed) return pointed;
            }
            return e && e.target ? e.target : null;
        }
        function _optionIndex(selectEl, opt){
            if (!selectEl || !opt || !selectEl.options) return -1;
            return Array.prototype.indexOf.call(selectEl.options, opt);
        }
        function _optionFromIndex(selectEl, idx){
            if (!selectEl || !selectEl.options || idx < 0 || idx >= selectEl.options.length) return null;
            return selectEl.options[idx];
        }
        function _optionVisualHeight(selectEl){
            if (!selectEl || !selectEl.options || !selectEl.options.length) return 0;
            for (var i = 0; i < selectEl.options.length; i++) {
                var rect = selectEl.options[i].getBoundingClientRect ? selectEl.options[i].getBoundingClientRect() : null;
                if (rect && rect.height > 0) return rect.height;
            }
            return 0;
        }
        function _selectedOption(selectEl){
            if (!selectEl || !selectEl.options || !selectEl.options.length) return null;
            if (selectEl.selectedOptions && selectEl.selectedOptions.length) {
                return selectEl.selectedOptions[0];
            }
            if (typeof selectEl.selectedIndex === 'number' && selectEl.selectedIndex >= 0) {
                return _optionFromIndex(selectEl, selectEl.selectedIndex);
            }
            return null;
        }
        function _optionIndexAtPoint(selectEl, e){
            if (!selectEl || !selectEl.options || !selectEl.options.length) return -1;
            var node = _eventPointNode(e);
            if (node && node.tagName === 'OPTION') return _optionIndex(selectEl, node);
            if (node && node.closest) {
                var optNode = node.closest('option');
                if (optNode) return _optionIndex(selectEl, optNode);
            }
            if (!e || typeof e.clientX !== 'number' || typeof e.clientY !== 'number') return -1;
            var rect = selectEl.getBoundingClientRect();
            if (
                e.clientX < rect.left || e.clientX > rect.right ||
                e.clientY < rect.top || e.clientY > rect.bottom
            ) {
                return -1;
            }
            var optionCount = selectEl.options.length;
            var optionHeight = _optionVisualHeight(selectEl);
            if ((!optionHeight || !isFinite(optionHeight)) && selectEl.scrollHeight && selectEl.scrollHeight > 0) {
                optionHeight = selectEl.scrollHeight / optionCount;
            }
            if ((!optionHeight || !isFinite(optionHeight)) && rect.height > 0) {
                var visibleRows = Math.max(1, Math.min(optionCount, Number(selectEl.size) || optionCount));
                optionHeight = rect.height / visibleRows;
            }
            if (!optionHeight || !isFinite(optionHeight)) return -1;
            var localY = (e.clientY - rect.top) + (selectEl.scrollTop || 0);
            var contentHeight = optionHeight * optionCount;
            if (localY < 0 || localY >= contentHeight) return -1;
            var rawIdx = Math.floor(localY / optionHeight);
            if (rawIdx < 0) rawIdx = 0;
            if (rawIdx >= optionCount) rawIdx = optionCount - 1;
            return rawIdx;
        }
        function _optionAtPoint(selectEl, e){
            return _optionFromIndex(selectEl, _optionIndexAtPoint(selectEl, e));
        }
        function _activeOptionForEvent(selectEl, e){
            return _optionAtPoint(selectEl, e) || _selectedOption(selectEl);
        }
        function _selectAtPoint(e){
            var node = _eventPointNode(e);
            if (!node || !node.closest) return null;
            return node.closest('.calc-file-list select');
        }
        function _batchSelectAtPoint(e){
            var node = _eventPointNode(e);
            if (!node || !node.closest) return null;
            return node.closest('.calc-xyz-batch-select select');
        }
        function _leftPaneAtPoint(e){
            var node = _eventPointNode(e);
            if (!node || !node.closest) return null;
            return node.closest('.calc-left');
        }
        function _pointInsideElement(el, e){
            if (!el || !e || typeof e.clientX !== 'number' || typeof e.clientY !== 'number') return false;
            var rect = el.getBoundingClientRect();
            return (
                e.clientX >= rect.left &&
                e.clientX <= rect.right &&
                e.clientY >= rect.top &&
                e.clientY <= rect.bottom
            );
        }
        function _hasExternalFiles(e){
            var dt = e && e.dataTransfer;
            if (!dt) return false;
            if (dt.files && dt.files.length) return true;
            var types = Array.prototype.slice.call(dt.types || []);
            return types.indexOf('Files') >= 0;
        }
        function _setDropActive(selectEl, active){
            if (!selectEl || !selectEl.classList) return;
            if (active) selectEl.classList.add('calc-drop-active');
            else selectEl.classList.remove('calc-drop-active');
        }
        function _widgetField(root, cls){
            if (!root) return null;
            return root.querySelector('.' + cls + ' textarea, .' + cls + ' input');
        }
        function _getWidgetFieldValue(root, cls){
            var field = _widgetField(root, cls);
            return field ? String(field.value == null ? '' : field.value) : '';
        }
        function _setWidgetField(root, cls, value){
            var field = _widgetField(root, cls);
            if (!field) return false;
            var strVal = String(value == null ? '' : value);
            var nativeSet = Object.getOwnPropertyDescriptor(
                Object.getPrototypeOf(field), 'value'
            ) || Object.getOwnPropertyDescriptor(HTMLInputElement.prototype, 'value')
              || Object.getOwnPropertyDescriptor(HTMLTextAreaElement.prototype, 'value');
            if (nativeSet && nativeSet.set) {
                nativeSet.set.call(field, strVal);
            } else {
                field.value = strVal;
            }
            field.dispatchEvent(new Event('input', { bubbles: true }));
            field.dispatchEvent(new Event('change', { bubbles: true }));
            return true;
        }
        function _setWidgetInput(root, cls, value){
            return _setWidgetField(root, cls, value);
        }
        function _clickWidgetButton(root, cls){
            if (!root) return false;
            var node = root.querySelector('.' + cls);
            if (!node) return false;
            var btn = null;
            if (node.tagName === 'BUTTON') {
                btn = node;
            } else if (node.querySelector) {
                btn = node.querySelector('button');
            }
            if (!btn && typeof node.click === 'function') {
                node.click();
                return true;
            }
            if (!btn) return false;
            btn.click();
            return true;
        }
        function _cookieValue(name){
            var parts = document.cookie ? document.cookie.split(';') : [];
            for (var i = 0; i < parts.length; i++) {
                var raw = String(parts[i] || '').trim();
                if (raw.indexOf(name + '=') === 0) {
                    return decodeURIComponent(raw.slice(name.length + 1));
                }
            }
            return '';
        }
        function _pageConfig(){
            try {
                var node = document.getElementById('jupyter-config-data');
                return node ? JSON.parse(node.textContent || '{}') : {};
            } catch (_err) {
                return {};
            }
        }
        function _baseUrl(){
            var cfg = _pageConfig();
            return String(cfg.baseUrl || '/');
        }
        function _uploadStagingRootRel(root){
            if (!root || !window.__delfinCalcUploadStagingRoots) return '';
            var classes = Array.prototype.slice.call(root.classList || []);
            for (var i = 0; i < classes.length; i++) {
                var cls = classes[i];
                if (cls.indexOf('calc-scope-') === 0 && window.__delfinCalcUploadStagingRoots[cls]) {
                    return String(window.__delfinCalcUploadStagingRoots[cls] || '');
                }
            }
            return '';
        }
        function _currentExplorerRelPath(root){
            if (!root) return '';
            var label = root.querySelector('.calc-path-label');
            var text = label ? String(label.textContent || label.innerText || '') : '';
            var idx = text.indexOf('Path:');
            if (idx >= 0) text = text.slice(idx + 5);
            text = text.replace(/\u00a0/g, ' ').trim();
            if (!text || text === '/') return '';
            return text.replace(/^\/+/, '').replace(/\/+$/, '');
        }
        function _setOpsStatus(root, message, color){
            if (!root) return;
            var node = root.querySelector('.calc-ops-status .widget-html-content') || root.querySelector('.calc-ops-status');
            if (!node) return;
            node.innerHTML = '<span style="color:' + String(color || '#555') + ';">' + String(message || '') + '</span>';
        }
        function _encodeContentsPath(path){
            return String(path || '')
                .split('/')
                .filter(function(part){ return !!part; })
                .map(function(part){ return encodeURIComponent(part); })
                .join('/');
        }
        function _normalizeRelParts(path){
            var raw = String(path || '').replace(/\\/g, '/');
            var parts = raw.split('/').filter(function(part){ return !!part && part !== '.'; });
            if (!parts.length) return [];
            for (var i = 0; i < parts.length; i++) {
                if (parts[i] === '..') throw new Error('Upload path cannot contain ".."');
            }
            return parts;
        }
        function _joinRelParts(parts){
            return parts.filter(function(part){ return !!part; }).join('/');
        }
        function _mergeRelParts(){
            var merged = [];
            for (var i = 0; i < arguments.length; i++) {
                var segment = arguments[i];
                if (!segment) continue;
                merged = merged.concat(_normalizeRelParts(segment));
            }
            return merged;
        }
        function _basename(path){
            var parts = _normalizeRelParts(path);
            return parts.length ? parts[parts.length - 1] : '';
        }
        function _widgetText(root, cls){
            if (!root) return '';
            var node = root.querySelector('.' + cls);
            return node ? String(node.textContent || node.innerText || '').trim() : '';
        }
        function _currentUploadAck(root){
            return _getWidgetFieldValue(root, 'calc-upload-ack') || _widgetText(root, 'calc-upload-ack-label');
        }
        function _sleep(ms){
            return new Promise(function(resolve){ setTimeout(resolve, ms); });
        }
        function _waitForUploadAck(root, previousAck, timeoutMs){
            return new Promise(function(resolve, reject){
                var start = Date.now();
                var prev = String(previousAck == null ? '' : previousAck);
                function check(){
                    var cur = _currentUploadAck(root);
                    if (cur && cur !== prev) {
                        resolve(cur);
                        return true;
                    }
                    if ((Date.now() - start) >= timeoutMs) {
                        reject(new Error('Explorer upload bridge timed out.'));
                        return true;
                    }
                    return false;
                }
                if (check()) return;
                var timer = setInterval(function(){
                    if (check()) clearInterval(timer);
                }, 40);
            });
        }
        async function _sendUploadChunk(root, meta, chunkData){
            var prevAck = _currentUploadAck(root);
            if (!_setWidgetInput(root, 'calc-upload-meta', JSON.stringify(meta))) {
                throw new Error('Upload meta widget not found');
            }
            if (!_setWidgetInput(root, 'calc-upload-chunk', chunkData)) {
                throw new Error('Upload chunk widget not found');
            }
            if (_clickWidgetButton(root, 'calc-upload-trigger-btn')) {
                await _waitForUploadAck(root, prevAck, 300000);
                await _sleep(5);
                return;
            }
            root._delfinUploadSeq = Math.max(
                Number(root._delfinUploadSeq) || 0,
                parseInt(String(_currentUploadAck(root) || '0'), 10) || 0
            ) + 1;
            if (!_setWidgetInput(root, 'calc-upload-seq', String(root._delfinUploadSeq))) {
                throw new Error('Upload trigger controls not found');
            }
            await _waitForUploadAck(root, prevAck, 300000);
            await _sleep(5);
        }
        async function _putContentsModel(relPath, model){
            var url = _baseUrl().replace(/\/?$/, '/') + 'api/contents/' + _encodeContentsPath(relPath);
            var headers = { 'Content-Type': 'application/json' };
            var xsrf = _cookieValue('_xsrf');
            if (xsrf) headers['X-XSRFToken'] = xsrf;
            var response = await fetch(url, {
                method: 'PUT',
                credentials: 'same-origin',
                headers: headers,
                body: JSON.stringify(model)
            });
            if (!response.ok) {
                var body = await response.text();
                throw new Error('Upload failed (' + response.status + '): ' + body.slice(0, 240));
            }
            return response;
        }
        async function _ensureContentsDirectory(relDir){
            var parts = _normalizeRelParts(relDir);
            var current = [];
            for (var i = 0; i < parts.length; i++) {
                current.push(parts[i]);
                await _putContentsModel(_joinRelParts(current), { type: 'directory' });
            }
        }
        function _hideCtxMenu(){
            var menu = document.getElementById('delfin-explorer-ctx-menu');
            if (menu) menu.style.display = 'none';
        }
        function _ensureCtxMenu(){
            var menu = document.getElementById('delfin-explorer-ctx-menu');
            if (menu) return menu;
            menu = document.createElement('div');
            menu.id = 'delfin-explorer-ctx-menu';
            menu.style.cssText = [
                'position:fixed',
                'z-index:2147483647',
                'display:none',
                'min-width:170px',
                'background:#fff',
                'border:1px solid #c8c8c8',
                'box-shadow:0 8px 24px rgba(0,0,0,0.18)',
                'border-radius:6px',
                'padding:4px'
            ].join(';');

            function makeItem(label, onClick){
                var item = document.createElement('button');
                item.type = 'button';
                item.textContent = label;
                item.style.cssText = [
                    'display:block',
                    'width:100%',
                    'text-align:left',
                    'padding:7px 10px',
                    'border:none',
                    'background:transparent',
                    'cursor:pointer',
                    'border-radius:4px',
                    'font:13px sans-serif'
                ].join(';');
                item.addEventListener('mouseenter', function(){ item.style.background = '#f2f6ff'; });
                item.addEventListener('mouseleave', function(){ item.style.background = 'transparent'; });
                item.addEventListener('click', function(ev){
                    ev.preventDefault();
                    ev.stopPropagation();
                    try { onClick(menu); } finally { _hideCtxMenu(); }
                });
                return item;
            }

            menu.appendChild(makeItem('New Folder', function(menuEl){
                var root = menuEl._root;
                _clickWidgetButton(root, 'calc-cmd-explorer-new-btn');
            }));
            menu.appendChild(makeItem('Rename', function(menuEl){
                var root = menuEl._root;
                var opt = menuEl._ctxOption;
                if (!root || !opt) return;
                var currentName = _labelToName(_labelText(opt));
                if (!currentName) return;
                var nextName = window.prompt('Rename to:', currentName);
                if (nextName === null) return;
                _setWidgetInput(root, 'calc-cmd-action-source', currentName);
                _setWidgetInput(root, 'calc-cmd-rename-input', nextName);
                _clickWidgetButton(root, 'calc-cmd-rename-btn');
            }));
            document.body.appendChild(menu);
            return menu;
        }
        function _triggerMove(root, sourceOpt, targetOpt){
            if (!root || !sourceOpt || !targetOpt || sourceOpt === targetOpt) return;
            var targetLabel = _labelText(targetOpt);
            if (!_isFolderLabel(targetLabel)) return;
            var targetName = _labelToName(targetLabel);
            if (!targetName) return;
            var sourceName = _labelToName(_labelText(sourceOpt));
            if (!sourceName) return;
            _setWidgetInput(root, 'calc-cmd-action-source', sourceName);
            _setWidgetInput(root, 'calc-cmd-move-target', targetName);
            setTimeout(function(){
                _clickWidgetButton(root, 'calc-cmd-move-btn');
            }, 20);
        }
        function _entryFile(entry){
            return new Promise(function(resolve, reject){
                try { entry.file(resolve, reject); }
                catch (err) { reject(err); }
            });
        }
        function _readDirEntries(reader){
            return new Promise(function(resolve, reject){
                try { reader.readEntries(resolve, reject); }
                catch (err) { reject(err); }
            });
        }
        async function _walkDroppedEntry(entry, prefix, files){
            if (!entry) return;
            if (entry.isFile) {
                var file = await _entryFile(entry);
                files.push({
                    file: file,
                    relativePath: String(prefix || '') + String(file.name || 'upload.bin'),
                });
                return;
            }
            if (!entry.isDirectory) return;
            var nextPrefix = String(prefix || '') + String(entry.name || 'folder') + '/';
            var reader = entry.createReader();
            while (true) {
                var entries = await _readDirEntries(reader);
                if (!entries || !entries.length) break;
                for (var i = 0; i < entries.length; i++) {
                    await _walkDroppedEntry(entries[i], nextPrefix, files);
                }
            }
        }
        async function _collectDroppedFiles(e){
            var dt = e && e.dataTransfer;
            if (!dt) return [];
            var files = [];
            var items = Array.prototype.slice.call(dt.items || []);
            var usedEntries = false;
            for (var i = 0; i < items.length; i++) {
                var item = items[i];
                if (!item || typeof item.webkitGetAsEntry !== 'function') continue;
                var entry = item.webkitGetAsEntry();
                if (!entry) continue;
                usedEntries = true;
                await _walkDroppedEntry(entry, '', files);
            }
            if (!usedEntries) {
                Array.prototype.forEach.call(dt.files || [], function(file){
                    if (!file) return;
                    files.push({
                        file: file,
                        relativePath: file.webkitRelativePath || file.name || 'upload.bin',
                    });
                });
            }
            return files;
        }
        function _readFileAsBase64(file){
            return new Promise(function(resolve, reject){
                var reader = new FileReader();
                reader.onload = function(){
                    var result = String(reader.result || '');
                    var comma = result.indexOf(',');
                    resolve(comma >= 0 ? result.slice(comma + 1) : result);
                };
                reader.onerror = function(){
                    reject(reader.error || new Error('Failed to read dropped file'));
                };
                reader.readAsDataURL(file);
            });
        }
        function _clipboardFiles(e){
            var cd = e && e.clipboardData;
            if (!cd) return [];
            var files = [];
            var items = Array.prototype.slice.call(cd.items || []);
            items.forEach(function(item){
                if (!item || item.kind !== 'file') return;
                var file = typeof item.getAsFile === 'function' ? item.getAsFile() : null;
                if (!file) return;
                files.push({
                    file: file,
                    relativePath: file.name || 'upload.bin',
                });
            });
            if (!files.length) {
                Array.prototype.forEach.call(cd.files || [], function(file){
                    if (!file) return;
                    files.push({
                        file: file,
                        relativePath: file.name || 'upload.bin',
                    });
                });
            }
            return files;
        }
        function _targetFolderName(selectEl, e, preferSelected){
            var targetOpt = e ? _optionAtPoint(selectEl, e) : null;
            if (!targetOpt && preferSelected && selectEl && selectEl.selectedOptions && selectEl.selectedOptions.length) {
                targetOpt = selectEl.selectedOptions[0];
            }
            var targetLabel = _labelText(targetOpt);
            return _isFolderLabel(targetLabel) ? _labelToName(targetLabel) : '';
        }
        async function _uploadBrowserFiles(root, browserFiles, targetName, batchPrefix){
            console.log('[DELFIN] _uploadBrowserFiles called, files:', browserFiles.length, 'targetName:', targetName);
            if (!root || !browserFiles || !browserFiles.length) return false;
            _setOpsStatus(root, 'Uploading ' + browserFiles.length + ' file(s)...', '#555');
            var chunkChars = 1000000;
            var batchId = String(batchPrefix || 'drop') + '_' + Date.now() + '_' + Math.random().toString(36).slice(2, 10);
            var target = String(targetName || '');
            for (var i = 0; i < browserFiles.length; i++) {
                var browserEntry = browserFiles[i];
                var fileObj = browserEntry && browserEntry.file ? browserEntry.file : null;
                if (!fileObj) continue;
                var relativePath = String(browserEntry.relativePath || fileObj.webkitRelativePath || fileObj.name || ('upload_' + i));
                var payload = await _readFileAsBase64(fileObj);
                var totalChunks = Math.max(1, Math.ceil(payload.length / chunkChars));
                for (var chunkIndex = 0; chunkIndex < totalChunks; chunkIndex++) {
                    var chunkData = payload.slice(chunkIndex * chunkChars, (chunkIndex + 1) * chunkChars);
                    await _sendUploadChunk(root, {
                        batch_id: batchId,
                        upload_id: batchId + ':' + i,
                        batch_total: browserFiles.length,
                        chunk_index: chunkIndex,
                        chunk_total: totalChunks,
                        target: target,
                        name: _basename(relativePath) || (fileObj.name || ('upload_' + i)),
                        relative_path: relativePath,
                    }, chunkData);
                }
            }
            _setOpsStatus(root, 'Upload sent. Processing...', '#2e7d32');
            return true;
        }
        async function _uploadDroppedFiles(root, selectEl, e){
            if (!root || !selectEl || !_hasExternalFiles(e)) return false;
            var droppedFiles = await _collectDroppedFiles(e);
            if (!droppedFiles.length) return false;
            return _uploadBrowserFiles(root, droppedFiles, _targetFolderName(selectEl, e, false), 'drop');
        }
        async function _uploadPastedFiles(root, selectEl, e){
            if (!root || !selectEl) return false;
            var pastedFiles = _clipboardFiles(e);
            if (!pastedFiles.length) return false;
            return _uploadBrowserFiles(root, pastedFiles, _targetFolderName(selectEl, null, true), 'paste');
        }
        function _installExplorerSelect(selectEl){
            if (!selectEl || selectEl._delfinExplorerReady) return;
            selectEl._delfinExplorerReady = true;
            var root = selectEl.closest('.calc-tab');
            if (!root) return;

            var dragState = { src: null, x: 0, y: 0, moved: false };
            var ctxOption = null;
            function _rememberExplorerScroll(){
                root._delfinExplorerScrollTop = selectEl.scrollTop || 0;
                root._delfinExplorerKeepScroll = true;
            }
            function _restoreExplorerScroll(){
                if (!root._delfinExplorerKeepScroll) return;
                var top = Number(root._delfinExplorerScrollTop) || 0;
                function apply(){
                    if (!selectEl || !selectEl.isConnected) return;
                    if (root.querySelector('.calc-file-list select') !== selectEl) return;
                    selectEl.scrollTop = top;
                }
                apply();
                requestAnimationFrame(apply);
                setTimeout(apply, 0);
                setTimeout(function(){
                    apply();
                    if (selectEl && selectEl.isConnected && root.querySelector('.calc-file-list select') === selectEl) {
                        root._delfinExplorerKeepScroll = false;
                    }
                }, 120);
            }
            function _armExplorerDblClick(opt, e){
                var label = _labelText(opt);
                if (!label || label.charAt(0) === '(') {
                    root._dblLastTime = 0;
                    root._dblLastLabel = '';
                    return;
                }
                root._dblLastTime = Date.now();
                root._dblLastLabel = label;
                root._dblLastX = (e && typeof e.clientX === 'number') ? e.clientX : 0;
                root._dblLastY = (e && typeof e.clientY === 'number') ? e.clientY : 0;
            }

            function refreshDraggable(){
                var opts = Array.prototype.slice.call(selectEl.options || []);
                opts.forEach(function(opt){
                    var txt = _labelText(opt);
                    opt.draggable = txt && txt.charAt(0) !== '(';
                });
            }
            _restoreExplorerScroll();
            selectEl.addEventListener('scroll', function(){
                root._delfinExplorerScrollTop = selectEl.scrollTop || 0;
            }, true);

            selectEl.addEventListener('contextmenu', function(e){
                var menu = _ensureCtxMenu();
                var opt = _optionAtPoint(selectEl, e);
                e.preventDefault();
                e.stopPropagation();
                if (opt) {
                    ctxOption = opt;
                } else {
                    ctxOption = null;
                }
                menu._root = root;
                menu._ctxOption = ctxOption;
                menu.style.left = (e.clientX + 2) + 'px';
                menu.style.top = (e.clientY + 2) + 'px';
                menu.style.display = 'block';
            }, true);

            /* --- Click selection + drag tracking --- */
            var lastIdx = -1;

            selectEl.addEventListener('mousedown', function(e){
                if (e.button != null && e.button !== 0) return;
                var opt = _optionAtPoint(selectEl, e);
                if (!opt) {
                    dragState.src = null;
                    dragState.moved = false;
                    root._dblLastTime = 0;
                    root._dblLastLabel = '';
                    root._dblLastX = 0;
                    root._dblLastY = 0;
                    _rememberExplorerScroll();
                    e.preventDefault();
                    Array.prototype.forEach.call(selectEl.options || [], function(item){
                        item.selected = false;
                    });
                    selectEl.dispatchEvent(new Event('change', { bubbles: true }));
                    selectEl.focus();
                    _restoreExplorerScroll();
                    return;
                }
                /* Track drag start */
                dragState.src = opt;
                dragState.x = e.clientX || 0;
                dragState.y = e.clientY || 0;
                dragState.moved = false;

                /* Prevent native select behavior – we manage selection ourselves */
                e.preventDefault();

                var opts = Array.from(selectEl.options);
                var idx = _optionIndex(selectEl, opt);
                if (idx < 0) {
                    dragState.src = null;
                    return;
                }

                /* Selection logic */
                if (e.shiftKey && lastIdx >= 0) {
                    var lo = Math.min(lastIdx, idx);
                    var hi = Math.max(lastIdx, idx);
                    if (!e.ctrlKey && !e.metaKey) {
                        for (var i = 0; i < opts.length; i++) opts[i].selected = false;
                    }
                    for (var j = lo; j <= hi; j++) opts[j].selected = true;
                } else if (e.ctrlKey || e.metaKey) {
                    opt.selected = !opt.selected;
                } else {
                    for (var k = 0; k < opts.length; k++) opts[k].selected = false;
                    opt.selected = true;
                }
                lastIdx = idx;
                _armExplorerDblClick(opt, e);
                _rememberExplorerScroll();

                selectEl.dispatchEvent(new Event('change', { bubbles: true }));
                selectEl.focus();
                _restoreExplorerScroll();
            }, true);

            selectEl.addEventListener('mousemove', function(e){
                if (!dragState.src) return;
                var dx = Math.abs((e.clientX || 0) - dragState.x);
                var dy = Math.abs((e.clientY || 0) - dragState.y);
                if (dx + dy > 8) dragState.moved = true;
            }, true);

            selectEl.addEventListener('mouseup', function(e){
                var targetOpt = _optionAtPoint(selectEl, e);
                if (!dragState.src || !dragState.moved || !targetOpt) {
                    dragState.src = null;
                    dragState.moved = false;
                    return;
                }
                e.preventDefault();
                _triggerMove(root, dragState.src, targetOpt);
                dragState.src = null;
                dragState.moved = false;
            }, true);

            /* --- Keyboard shortcuts --- */
            selectEl.addEventListener('keydown', function(e){
                var action = '';
                if ((e.ctrlKey || e.metaKey) && e.key === 'c') action = 'copy';
                else if ((e.ctrlKey || e.metaKey) && e.key === 'x') action = 'cut';
                else if ((e.ctrlKey || e.metaKey) && e.key === 'v') action = 'paste';
                else if (e.key === 'Delete' || e.key === 'Backspace') action = 'delete';
                else if (e.key === 'Enter') action = 'open';
                else if (e.key === 'Escape') action = 'deselect';
                else if (e.key === 'F2') action = 'rename';
                if (action) {
                    e.preventDefault();
                    e.stopPropagation();
                    _setWidgetInput(root, 'calc-cmd-keyboard-action', action);
                }
            }, true);

            selectEl.addEventListener('dragstart', function(e){
                var opt = _activeOptionForEvent(selectEl, e);
                if (!opt) return;
                selectEl._delfinDragSource = opt;
                try {
                    e.dataTransfer.effectAllowed = 'move';
                    e.dataTransfer.setData('text/plain', _labelText(opt));
                } catch (_err) {}
            }, true);

            selectEl.addEventListener('dragenter', function(e){
                if (!_hasExternalFiles(e)) return;
                e.preventDefault();
                _setDropActive(selectEl, true);
            }, true);

            selectEl.addEventListener('dragover', function(e){
                if (_hasExternalFiles(e)) {
                    e.preventDefault();
                    _setDropActive(selectEl, true);
                    try { e.dataTransfer.dropEffect = 'copy'; } catch (_err) {}
                    return;
                }
                var opt = _optionAtPoint(selectEl, e);
                if (opt && _isFolderLabel(_labelText(opt))) e.preventDefault();
            }, true);

            selectEl.addEventListener('dragleave', function(e){
                if (!_hasExternalFiles(e)) return;
                setTimeout(function(){ _setDropActive(selectEl, false); }, 0);
            }, true);

            selectEl.addEventListener('drop', function(e){
                if (_hasExternalFiles(e)) {
                    e.preventDefault();
                    e.stopPropagation();
                    _setDropActive(selectEl, false);
                    selectEl._delfinDragSource = null;
                    _uploadDroppedFiles(root, selectEl, e).catch(function(err){
                        _setOpsStatus(root, 'Explorer upload failed: ' + String(err && err.message ? err.message : err), '#d32f2f');
                        console.error('Explorer upload failed', err);
                    });
                    return;
                }
                var targetOpt = _optionAtPoint(selectEl, e);
                var sourceOpt = selectEl._delfinDragSource || null;
                _setDropActive(selectEl, false);
                if (!sourceOpt || !targetOpt) {
                    selectEl._delfinDragSource = null;
                    return;
                }
                e.preventDefault();
                _triggerMove(root, sourceOpt, targetOpt);
                selectEl._delfinDragSource = null;
            }, true);

            refreshDraggable();
            try {
                new MutationObserver(refreshDraggable).observe(selectEl, {
                    childList: true,
                    subtree: true,
                });
            } catch (_err) {}
        }
        function _installBatchSelect(selectEl){
            if (!selectEl || selectEl._delfinBatchSelectReady) return;
            selectEl._delfinBatchSelectReady = true;
            var root = selectEl.closest('.calc-tab');
            if (!root) return;

            function _isBatchToggleZone(e){
                if (!e || typeof e.clientX !== 'number') return false;
                var opt = _optionAtPoint(selectEl, e);
                if (!opt || !opt.getBoundingClientRect) return false;
                var rect = opt.getBoundingClientRect();
                var localX = e.clientX - rect.left;
                return localX >= 0 && localX <= 64;
            }
            function _styleBatchOption(opt){
                if (!opt) return;
                var text = _labelText(opt);
                if (text.indexOf('✔ ') === 0) {
                    opt.style.color = '#2e7d32';
                    opt.style.fontWeight = '700';
                } else {
                    opt.style.color = '';
                    opt.style.fontWeight = '';
                }
            }
            function _applyBatchOptionStyles(){
                Array.prototype.forEach.call(selectEl.options || [], function(opt){
                    _styleBatchOption(opt);
                });
            }
            function _toggleBatchOptionLabel(opt){
                if (!opt) return;
                var text = _labelText(opt);
                if (!text) return;
                if (text.indexOf('□ ') === 0) {
                    opt.textContent = text.replace(/^□\s+/, '✔ ');
                    _styleBatchOption(opt);
                    return;
                }
                if (text.indexOf('✔ ') === 0) {
                    opt.textContent = text.replace(/^✔\s+/, '□ ');
                    _styleBatchOption(opt);
                }
            }

            function _resetBatchDblClick(){
                root._batchDblLastTime = 0;
                root._batchDblLastLabel = '';
                root._batchDblLastX = 0;
                root._batchDblLastY = 0;
            }
            function _armBatchDblClick(opt, e){
                var label = _labelText(opt);
                if (!label || label.charAt(0) === '(') {
                    _resetBatchDblClick();
                    return;
                }
                root._batchDblLastTime = Date.now();
                root._batchDblLastLabel = label;
                root._batchDblLastX = (e && typeof e.clientX === 'number') ? e.clientX : 0;
                root._batchDblLastY = (e && typeof e.clientY === 'number') ? e.clientY : 0;
            }

            selectEl.addEventListener('mousedown', function(e){
                if (e.button != null && e.button !== 0) return;
                var opt = _optionAtPoint(selectEl, e);
                if (!opt) {
                    _resetBatchDblClick();
                    return;
                }
                if (_isBatchToggleZone(e)) {
                    _resetBatchDblClick();
                    e.preventDefault();
                    e.stopPropagation();
                    _toggleBatchOptionLabel(opt);
                    _setWidgetInput(root, 'calc-cmd-xyz-batch-toggle', 'idx:' + String(_optionIndex(selectEl, opt)));
                    return;
                }
                _armBatchDblClick(opt, e);
            }, true);

            _applyBatchOptionStyles();
            try {
                new MutationObserver(_applyBatchOptionStyles).observe(selectEl, {
                    childList: true,
                    subtree: true,
                });
            } catch (_err) {}
        }
        function _installExternalFileDrop(root){
            if (!root || root._delfinExternalDropReady) return;
            console.log('[DELFIN] _installExternalFileDrop on', root.className);
            root._delfinExternalDropReady = true;
            root._delfinPasteArmed = false;

            /* --- Double-click via delegated mousedown on stable root ---
             * ipywidgets may re-create the <select> between clicks, so native
             * dblclick never fires.  We therefore remember the first click on
             * the stable root and only use the second click to confirm timing
             * and pointer position. */
            root._dblLastTime = 0;
            root._dblLastLabel = '';
            root._dblLastX = 0;
            root._dblLastY = 0;
            root._batchDblLastTime = 0;
            root._batchDblLastLabel = '';
            root._batchDblLastX = 0;
            root._batchDblLastY = 0;
            root.addEventListener('mousedown', function(e){
                if (e.button != null && e.button !== 0) return;
                var sel = _selectAtPoint(e);
                if (sel && sel.closest('.calc-file-list')) {
                    var currentOpt = _optionAtPoint(sel, e);
                    if (!currentOpt) return;
                    var currentLabel = _labelText(currentOpt);
                    var label = String(root._dblLastLabel || '');
                    if (!label || label.charAt(0) === '(') return;
                    var now = Date.now();
                    var dx = Math.abs((typeof e.clientX === 'number' ? e.clientX : 0) - (Number(root._dblLastX) || 0));
                    var dy = Math.abs((typeof e.clientY === 'number' ? e.clientY : 0) - (Number(root._dblLastY) || 0));
                    if (currentLabel === label && (now - root._dblLastTime) < 500 && dx <= 20 && dy <= 20) {
                        root._dblLastTime = 0;
                        root._dblLastLabel = '';
                        root._dblLastX = 0;
                        root._dblLastY = 0;
                        e.preventDefault();
                        e.stopPropagation();
                        _setWidgetInput(root, 'calc-cmd-dblclick', label);
                    }
                    return;
                }
                sel = _batchSelectAtPoint(e);
                if (!sel || !sel.closest('.calc-xyz-batch-select')) return;
                if (typeof e.clientX === 'number') {
                    var batchOptForZone = _optionAtPoint(sel, e);
                    if (batchOptForZone && batchOptForZone.getBoundingClientRect) {
                        var batchRect = batchOptForZone.getBoundingClientRect();
                        var batchLocalX = e.clientX - batchRect.left;
                        if (batchLocalX >= 0 && batchLocalX <= 64) return;
                    }
                }
                var batchOpt = _optionAtPoint(sel, e);
                if (!batchOpt) return;
                var batchLabel = _labelText(batchOpt);
                var prevBatchLabel = String(root._batchDblLastLabel || '');
                if (!prevBatchLabel || prevBatchLabel.charAt(0) === '(') return;
                var batchNow = Date.now();
                var batchDx = Math.abs((typeof e.clientX === 'number' ? e.clientX : 0) - (Number(root._batchDblLastX) || 0));
                var batchDy = Math.abs((typeof e.clientY === 'number' ? e.clientY : 0) - (Number(root._batchDblLastY) || 0));
                if (batchLabel === prevBatchLabel && (batchNow - root._batchDblLastTime) < 500 && batchDx <= 20 && batchDy <= 20) {
                    root._batchDblLastTime = 0;
                    root._batchDblLastLabel = '';
                    root._batchDblLastX = 0;
                    root._batchDblLastY = 0;
                    e.preventDefault();
                    e.stopPropagation();
                    _setWidgetInput(root, 'calc-cmd-xyz-batch-dblclick', prevBatchLabel);
                }
            }, true);

            function activeSelectFromEvent(e){
                if (!pointerInsideExplorer(e)) return null;
                var selectEl = _selectAtPoint(e);
                if (selectEl && root.contains(selectEl)) return selectEl;
                return root.querySelector('.calc-file-list select');
            }
            function activeSelectForPaste(){
                return root.querySelector('.calc-file-list select');
            }
            function pointerInsideExplorer(e){
                var leftPane = root.querySelector('.calc-left');
                if (!leftPane) return false;
                var leftPaneAtPoint = _leftPaneAtPoint(e);
                if (leftPaneAtPoint) return root.contains(leftPaneAtPoint);
                return _pointInsideElement(leftPane, e);
            }
            function clearHighlight(){
                var active = root.querySelector('.calc-file-list select.calc-drop-active');
                if (active) _setDropActive(active, false);
            }
            function handleExternalDrag(e){
                if (!_hasExternalFiles(e)) return false;
                if (!pointerInsideExplorer(e)) return false;
                var selectEl = activeSelectFromEvent(e);
                e.preventDefault();
                e.stopPropagation();
                if (selectEl) _setDropActive(selectEl, true);
                return true;
            }

            document.addEventListener('dragenter', function(e){
                if (_hasExternalFiles(e)) {
                    e.preventDefault();
                    e.stopPropagation();
                }
                handleExternalDrag(e);
            }, true);
            document.addEventListener('dragover', function(e){
                if (_hasExternalFiles(e)) {
                    e.preventDefault();
                    e.stopPropagation();
                    try { e.dataTransfer.dropEffect = 'copy'; } catch (_err) {}
                }
                if (!handleExternalDrag(e)) clearHighlight();
            }, true);
            document.addEventListener('drop', function(e){
                console.log('[DELFIN] drop event', '_hasExternalFiles=', _hasExternalFiles(e), 'types=', e.dataTransfer ? Array.from(e.dataTransfer.types) : 'none');
                if (!_hasExternalFiles(e)) return;
                if (!pointerInsideExplorer(e)) return;
                e.preventDefault();
                if (typeof e.stopImmediatePropagation === 'function') e.stopImmediatePropagation();
                else e.stopPropagation();
                clearHighlight();
                var selectEl = activeSelectFromEvent(e);
                console.log('[DELFIN] drop selectEl found:', !!selectEl, 'root:', !!root);
                if (selectEl) {
                    _uploadDroppedFiles(root, selectEl, e).catch(function(err){
                        _setOpsStatus(root, 'Explorer upload failed: ' + String(err && err.message ? err.message : err), '#d32f2f');
                        console.error('[DELFIN] Explorer upload failed', err);
                    });
                } else {
                    console.warn('[DELFIN] No selectEl found for drop upload');
                }
            }, true);
            root.addEventListener('mousedown', function(){
                root._delfinPasteArmed = true;
            }, true);
            root.addEventListener('focusin', function(){
                root._delfinPasteArmed = true;
            }, true);
            document.addEventListener('mousedown', function(e){
                if (!root.contains(e.target)) root._delfinPasteArmed = false;
            }, true);
            document.addEventListener('paste', function(e){
                var selectEl = activeSelectForPaste();
                if (!selectEl) return;
                var files = _clipboardFiles(e);
                if (!files.length) return;
                var activeInside = root._delfinPasteArmed || root.contains(document.activeElement);
                if (!activeInside) return;
                e.preventDefault();
                e.stopPropagation();
                _uploadPastedFiles(root, selectEl, e).catch(function(err){
                    _setOpsStatus(root, 'Explorer paste upload failed: ' + String(err && err.message ? err.message : err), '#d32f2f');
                    console.error('Explorer paste upload failed', err);
                });
            }, true);
            document.addEventListener('dragleave', function(_e){
                setTimeout(clearHighlight, 0);
            }, true);
            document.addEventListener('dragend', clearHighlight, true);
        }
        function _scanAndInstall(root){
            if (!root || !root.querySelectorAll) return;
            var selects = [];
            if (root.matches && root.matches('.calc-file-list select')) selects.push(root);
            Array.prototype.forEach.call(root.querySelectorAll('.calc-file-list select'), function(sel){
                selects.push(sel);
            });
            Array.prototype.forEach.call(selects, _installExplorerSelect);

            var batchSelects = [];
            if (root.matches && root.matches('.calc-xyz-batch-select select')) batchSelects.push(root);
            Array.prototype.forEach.call(root.querySelectorAll('.calc-xyz-batch-select select'), function(sel){
                batchSelects.push(sel);
            });
            Array.prototype.forEach.call(batchSelects, _installBatchSelect);

            var tabs = [];
            if (root.matches && root.matches('.calc-tab')) tabs.push(root);
            Array.prototype.forEach.call(root.querySelectorAll('.calc-tab'), function(tab){
                tabs.push(tab);
            });
            Array.prototype.forEach.call(tabs, _installExternalFileDrop);
        }

        _scanAndInstall(document.body);
        var _scanRetryCount = 0;
        function _retryScan(){
            _scanAndInstall(document.body);
            _scanRetryCount++;
            if (_scanRetryCount < 30) {
                setTimeout(_retryScan, 500);
            }
        }
        setTimeout(_retryScan, 300);
        new MutationObserver(function(mutations){
            mutations.forEach(function(m){
                Array.prototype.forEach.call(m.addedNodes || [], function(node){
                    if (node && node.nodeType === 1) _scanAndInstall(node);
                });
            });
            _scanAndInstall(document.body);
        }).observe(document.body, { childList: true, subtree: true });

        document.addEventListener('click', _hideCtxMenu, true);
        document.addEventListener('scroll', _hideCtxMenu, true);
        document.addEventListener('keydown', function(e){
            if (e.key === 'Escape') _hideCtxMenu();
        }, true);
    })();
    """
    # -- layout -------------------------------------------------------------
    _archive_nav_children = [calc_table_btn] if _is_archive_tab else []
    _archive_selection_children = []
    if not _is_archive_tab:
        if _remote_archive_enabled:
            _archive_selection_children.append(calc_ssh_transfer_btn)
        _archive_selection_children.append(calc_move_archive_btn)
    else:
        if _remote_archive_enabled:
            _archive_selection_children.append(calc_ssh_transfer_btn)
        _archive_selection_children.append(calc_back_to_calculations_btn)
    calc_nav_controls_row = widgets.HBox(
        [calc_back_btn, calc_home_btn, calc_refresh_btn, calc_delete_btn, *_archive_nav_children],
        layout=widgets.Layout(
            width='100%', overflow_x='hidden',
            justify_content='flex-start', gap='6px',
        ),
    )
    calc_nav_selection_row = widgets.HBox(
        [calc_explorer_new_btn, calc_explorer_rename_btn, calc_duplicate_btn, *_archive_selection_children],
        layout=widgets.Layout(
            width='100%', overflow_x='hidden',
            justify_content='flex-start', gap='6px',
            display='flex', flex_flow='row wrap',
        ),
    )
    calc_nav_bar = widgets.VBox([
        calc_path_label,
        calc_nav_controls_row,
        calc_nav_selection_row,
        calc_new_folder_prompt_row,
        calc_rename_prompt_row,
        calc_duplicate_prompt_row,
        calc_hidden_upload,
        widgets.VBox(
            [calc_new_folder_input, calc_new_folder_btn, calc_rename_input, calc_rename_btn,
             calc_move_target_input, calc_upload_target_input,
             calc_action_source_input, calc_move_btn,
             calc_upload_meta_input, calc_upload_chunk_input,
             calc_upload_seq_input, calc_upload_ack_input,
             calc_upload_trigger_btn, calc_upload_ack_label,
             calc_dblclick_input, calc_xyz_batch_dblclick_input, calc_xyz_batch_toggle_input,
             calc_keyboard_action_input],
            layout=widgets.Layout(display='none'),
        ),
    ], layout=widgets.Layout(width='100%', overflow_x='hidden'))

    calc_left = widgets.VBox([
        calc_nav_bar,
        calc_filter_sort_row,
        calc_file_list,
        calc_transfer_jobs_panel,
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
        '.calc-tab .calc-table-panel { flex:1 1 0 !important; min-height:0 !important;'
        ' overflow:hidden !important; }'
        '.calc-table-output { flex:1 1 0 !important; min-height:0 !important;'
        ' overflow-y:auto !important; overflow-x:auto !important; max-height:none !important; }'
        '.calc-tab .calc-table-panel .widget-text { overflow:visible !important; }'
        '.calc-tab .calc-table-panel .widget-dropdown, .calc-tab .calc-table-panel .widget-text { flex:0 0 auto !important; }'
        '.calc-tab .calc-table-panel .widget-dropdown, .calc-tab .calc-table-panel .widget-dropdown select { overflow:hidden !important; }'
        '.calc-tab .calc-table-top-row > .widget-text { flex:1 1 0 !important; min-width:0 !important; width:1px !important; }'
        '.calc-tab .calc-table-top-row > .widget-text input { width:100% !important; min-width:0 !important; }'
        '.calc-tab .calc-table-col-row > .widget-text { flex:1 1 0 !important; min-width:0 !important; width:1px !important; height:26px !important; min-height:26px !important; }'
        '.calc-tab .calc-table-col-row > .widget-text input { width:100% !important; min-width:0 !important; height:26px !important; line-height:26px !important; padding:0 6px !important; box-sizing:border-box !important; }'
        '.calc-tab .calc-table-panel .widget-text { height:26px !important; min-height:26px !important; }'
        '.calc-tab .calc-table-panel .widget-text input { width:100% !important; overflow-x:hidden !important; overflow-y:hidden !important; height:26px !important; line-height:26px !important; padding:0 6px !important; box-sizing:border-box !important; }'
        '.calc-tab .calc-table-panel .widget-dropdown { height:26px !important; min-height:26px !important; }'
        '.calc-tab .calc-table-panel .widget-dropdown select { height:26px !important; max-height:26px !important; line-height:26px !important; overflow:hidden !important; padding:0 6px !important; box-sizing:border-box !important; }'
        '.calc-tab .calc-table-panel input::-webkit-scrollbar { width:0; height:0; display:none; }'
        '.calc-tab .calc-table-panel input { scrollbar-width:none; overflow:hidden !important; height:26px !important; line-height:26px !important; padding:0 6px !important; box-sizing:border-box !important; }'
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
        '.calc-hidden-upload { position:absolute !important; width:1px !important; height:1px !important;'
        ' opacity:0 !important; overflow:hidden !important; pointer-events:none !important; }'
        '.calc-hidden-upload input, .calc-hidden-upload button'
        ' { width:1px !important; height:1px !important; opacity:0 !important; pointer-events:none !important; }'
        '.calc-left .calc-file-list select.calc-drop-active'
        ' { border:2px dashed #1976d2 !important; background:#edf4ff !important; }'
        '.calc-left .widget-select select'
        ' { text-overflow: ellipsis; white-space: nowrap; }'
        '.calc-path-label { gap:2px !important; align-items:stretch !important; }'
        '.calc-path-label > .widget-html { width:100% !important; margin:0 !important; }'
        '.calc-path-label > .widget-text'
        ' { flex:1 1 0 !important; min-width:0 !important; width:1px !important;'
        ' max-width:100% !important; margin:0 !important; }'
        '.calc-path-label > .widget-text input'
        ' { width:100% !important; height:24px !important; line-height:24px !important;'
        ' padding:2px 5px !important; margin:0 !important; border:0 !important;'
        ' border-radius:0 !important; box-shadow:none !important; background:#eee !important;'
        ' font-family:monospace !important; overflow:hidden !important;'
        ' text-overflow:ellipsis !important; white-space:nowrap !important; }'
        '.calc-left .widget-text { flex:1 1 auto !important; min-width:0 !important; width:auto !important; }'
        '.calc-left .widget-text input { width:100% !important; }'
        '.calc-filter-row { width:100% !important; display:flex !important; }'
        '.calc-filter-row > .widget-text { flex:1 1 auto !important; min-width:0 !important; width:auto !important; }'
        '.calc-filter-row > .widget-text input { width:100% !important; min-width:0 !important; }'
        '.calc-filter-row > .widget-dropdown { flex:0 0 90px !important; margin-left:auto !important; }'
        '.calc-left.calc-transfer-jobs-mode .calc-filter-row,'
        ' .calc-left.calc-transfer-jobs-mode .calc-file-list'
        ' { display:none !important; }'
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
        '.calc-png-btn button { display:inline-flex !important; align-items:center !important;'
        ' justify-content:center !important; text-align:center !important; padding:0 !important;'
        ' line-height:30px !important; }'
        '.calc-png-btn button i, .calc-png-btn .fa { display:none !important; }'
        '.calc-png-btn button > span, .calc-png-btn button .widget-button-label {'
        ' display:flex !important; align-items:center !important; justify-content:center !important;'
        ' width:100% !important; height:100% !important; margin:0 auto !important;'
        ' text-align:center !important; line-height:30px !important; padding:0 !important; }'
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

    calc_right_children = [
        widgets.HBox([
            calc_file_info,
            widgets.HBox(
                [
                    calc_transfer_jobs_btn,
                    calc_copy_path_btn,
                    calc_copy_btn,
                    calc_download_btn,
                    calc_report_btn,
                    calc_view_toggle,
                ],
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
        calc_xyz_batch_panel,
        calc_print_mode_panel,
        calc_mo_plot_panel,
        calc_delete_confirm,
        calc_delete_status,
        calc_mol_container,
        calc_content_label,
        calc_preselect_container,
        calc_move_archive_confirm,
        calc_move_archive_status,
        calc_ops_status,
    ]
    if _is_archive_tab:
        calc_right_children.append(calc_table_panel)
    calc_right_children.extend([
        calc_recalc_toolbar,
        calc_xyz_workflow_toolbar,
        calc_nmr_panel,
        calc_censo_nmr_panel,
        calc_content_toolbar,
        calc_chunk_hidden_row,
        calc_override_status,
        calc_content_area,
        calc_edit_area,
    ])
    calc_right = widgets.VBox(calc_right_children, layout=widgets.Layout(
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
    calc_nav_controls_row.add_class('calc-nav-controls-row')
    calc_nav_selection_row.add_class('calc-nav-selection-row')
    calc_rename_prompt_row.add_class('calc-rename-prompt-row')
    calc_duplicate_prompt_row.add_class('calc-duplicate-prompt-row')
    calc_content_area.add_class('calc-content-area')
    calc_table_panel.add_class('calc-table-panel')
    calc_table_output.add_class('calc-table-output')

    # Enable draggable splitter + dynamic mol-viewer resize + $3Dmol patch
    # Stored as a plain string; the CALLER (create_dashboard in __init__.py)
    # runs ALL tab init scripts in one ctx.run_js() call so that no tab's
    # clear_output() wipes another tab's init JS.
    _init_js = (
        _explorer_interactions_js
        + "\n"
        + f"""
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
            window.__delfinCalcUploadStagingRoots = window.__delfinCalcUploadStagingRoots || {{}};
            window.__delfinCalcUploadStagingRoots[scopeKey] = {json.dumps(CALC_BROWSER_UPLOAD_STAGING_SCOPE_REL)};

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
            if (typeof $3Dmol === 'undefined' || $3Dmol._calcResizePatched) return;
            var orig = $3Dmol.createViewer;
            var wrapper = function() {{
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
            try {{
                $3Dmol.createViewer = wrapper;
                $3Dmol._calcResizePatched = true;
            }} catch (_e1) {{
                try {{
                    Object.defineProperty($3Dmol, 'createViewer', {{
                        value: wrapper, writable: true, configurable: true
                    }});
                    $3Dmol._calcResizePatched = true;
                }} catch (_e2) {{
                    $3Dmol._calcResizePatched = true;
                }}
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
            var rightPanel = scopeRoot.querySelector('.calc-right');
            var mv = scopeRoot.querySelector('.calc-mol-viewer');
            if (!rightPanel || !mv || mv.offsetParent === null) return;
            var scopeKey = {json.dumps(calc_scope_id)};
            var savedView = null;
            var savedViewScope = null;
            window._calcMolViewStateByScope = window._calcMolViewStateByScope || {{}};
            window._calcMolViewScopeKeyByScope = window._calcMolViewScopeKeyByScope || {{}};
            var currentViewer =
                (window._calcMolViewerByScope && window._calcMolViewerByScope[scopeKey])
                || (window._calcTrajViewerByScope && window._calcTrajViewerByScope[scopeKey])
                || null;
            if (currentViewer && typeof currentViewer.getView === 'function') {{
                try {{
                    savedView = currentViewer.getView();
                    savedViewScope = window._calcMolViewScopeKeyByScope[scopeKey] || null;
                    if (savedView && savedViewScope) {{
                        window._calcMolViewStateByScope[savedViewScope] = savedView;
                    }}
                }} catch (_e) {{}}
            }}
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
                    if (!style || style.display === 'none' || style.visibility === 'hidden') {{
                        continue;
                    }}
                    var childRect = child.getBoundingClientRect();
                    if (childRect.height > 0) reservedBelow += childRect.height;
                }}
            }}
            var availH = rightRect.bottom - mvRect.top - reservedBelow - 10;
            var row = mv.closest('.calc-mol-view-row');
            var rowRect = row ? row.getBoundingClientRect() : rightRect;
            var tray = scopeRoot ? scopeRoot.querySelector('.calc-xyz-tray-controls') : null;
            var trayStyle = tray ? window.getComputedStyle(tray) : null;
            var trayVisible = !!(tray && trayStyle && trayStyle.display !== 'none');
            var trayWidth = trayVisible ? tray.getBoundingClientRect().width : 0;
            var availW = Math.max(120, rowRect.width - trayWidth - 16);
            var h = Math.floor(availH * {CALC_MOL_DYNAMIC_SCALE});
            var w = Math.floor(Math.min(h * 1.2, availW));
            if (h < 80) h = 80;
            if (w < 120) w = 120;
            mv.style.width = w + 'px';
            mv.style.height = h + 'px';
            var stage = mv.querySelector(
                '[id^="mol3d_"], [id^="calc_trj_viewer_"], [id^="3dmolviewer_"], [id^="calc_trj"]'
            );
            if (stage) {{
                stage.style.width = w + 'px';
                stage.style.height = h + 'px';
            }}
            var v = currentViewer;
            if (v && typeof v.resize === 'function') {{ v.resize(); v.render(); }}
            if (savedView && v && typeof v.setView === 'function') {{
                try {{
                    v.setView(savedView);
                    v.render();
                }} catch (_e) {{}}
            }}
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
    )

    def calc_set_root(root_dir):
        """Switch browser root directory and reset to top level."""
        try:
            new_root = Path(root_dir).expanduser()
        except Exception:
            new_root = Path(root_dir)
        new_root.mkdir(parents=True, exist_ok=True)
        ctx.calc_dir = new_root
        if _is_archive_tab:
            ctx.archive_dir = new_root
        state['current_path'] = ''
        calc_list_directory()

    def calc_set_primary_root(root_dir):
        """Update the Calculations root used by Archive move-back actions."""
        try:
            new_root = Path(root_dir).expanduser()
        except Exception:
            new_root = Path(root_dir)
        new_root.mkdir(parents=True, exist_ok=True)
        ctx.primary_calc_dir = new_root

    if _remote_archive_enabled:
        _calc_update_transfer_jobs_visibility()

    return tab_widget, {
        'calc_list_directory': calc_list_directory,
        'calc_set_root': calc_set_root,
        'calc_set_primary_root': calc_set_primary_root,
        'init_js': _init_js,
        # --- Agent-accessible widgets ---
        # Navigation
        'calc_path_input': calc_path_input,
        'calc_sort_dropdown': calc_sort_dropdown,
        'calc_folder_search': calc_folder_search,
        'calc_search_input': calc_search_input,
        # File operations
        'calc_new_folder_btn': calc_new_folder_btn,
        'calc_new_folder_input': calc_new_folder_input,
        'calc_rename_btn': calc_rename_btn,
        'calc_rename_input': calc_rename_input,
        'calc_duplicate_btn': calc_duplicate_btn,
        'calc_copy_btn': calc_copy_btn,
        'calc_copy_path_btn': calc_copy_path_btn,
        # Transfer / move
        'calc_move_archive_btn': calc_move_archive_btn,
        'calc_back_to_calculations_btn': calc_back_to_calculations_btn,
        'calc_ssh_transfer_btn': calc_ssh_transfer_btn,
        # Visualization & report
        'calc_view_toggle': calc_view_toggle,
        'calc_view_png_btn': calc_view_png_btn,
        'calc_xyz_png_btn': calc_xyz_png_btn,
        'calc_report_btn': calc_report_btn,
        # Extract Table
        'calc_table_btn': calc_table_btn,
        'calc_table_file_input': calc_table_file_input,
        'calc_table_scope_dd': calc_table_scope_dd,
        'calc_table_recursive_cb': calc_table_recursive_cb,
        'calc_table_decimal_comma_btn': calc_table_decimal_comma_btn,
        'calc_table_preset_name': calc_table_preset_name,
        'calc_table_preset_save_btn': calc_table_preset_save_btn,
        'calc_table_add_col_btn': calc_table_add_col_btn,
        'calc_table_run_btn': calc_table_run_btn,
        'calc_table_csv_btn': calc_table_csv_btn,
        'calc_table_output': calc_table_output,
        # Recalc & editor
        'calc_recalc_btn': calc_recalc_btn,
        'calc_recalc_time': calc_recalc_time,
        'calc_submit_recalc_btn': calc_submit_recalc_btn,
        'calc_edit_area': calc_edit_area,
        # Options dropdown (Recalc / Smart Recalc / Override / etc.)
        'calc_options_dropdown': calc_options_dropdown,
        'calc_override_input': calc_override_input,
        'calc_override_time': calc_override_time,
        'calc_override_btn': calc_override_btn,
        # File browser selection
        'calc_file_list': calc_file_list,
        'calc_update_options_dropdown': calc_update_options_dropdown,
        # Delete (exposed for blocking only — agent must NOT click)
        'calc_delete_btn': calc_delete_btn,
        # Batch from XYZ
        'xyz_batch_prepare_export': _calc_xyz_batch_prepare_export,
        'xyz_batch_state': state,
        'xyz_batch_select': calc_xyz_batch_select,
        'xyz_batch_filename': calc_xyz_batch_filename,
        'xyz_batch_refresh': _calc_refresh_xyz_batch_selector,
        'xyz_batch_show_panel': _calc_show_xyz_batch_panel,
    }
